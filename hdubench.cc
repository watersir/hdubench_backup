/*
	Discovering Hard Disk Physical Geometry through Microbenchmarking
	http://blog.stuffedcow.net/2019/09/hard-disk-geometry-microbenchmarking/
	Henry Wong, 2019/09/06

	2022/10/12: Fixed overflow in find_next_track_boundary. Add more sanity checks to RPM measurement.
*/

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/ioctl.h>
#include <fcntl.h>
#include <unistd.h>
#include <malloc.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <sys/time.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <getopt.h>
#include <vector>
#include <stdint.h>

#include <linux/fs.h>


// Copy-pasting Tiny Mersenne Twister random number generator from tinymt64.h/c
// http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/TINYMT/index.html


/*
 * tinymt64 internal state vector and parameters
 */
struct TINYMT64_T {
    uint64_t status[2];
    uint32_t mat1;
    uint32_t mat2;
    uint64_t tmat;
};

typedef struct TINYMT64_T tinymt64_t;

#ifndef UINT64_C
#define UINT64_C(X) (X ## ULL)
#endif

#define TINYMT64_SH0 12
#define TINYMT64_SH1 11
#define TINYMT64_SH8 8
#define TINYMT64_MASK UINT64_C(0x7fffffffffffffff)
#define MIN_LOOP 8

/**
 * This function changes internal state of tinymt64.
 * Users should not call this function directly.
 * @param random tinymt internal status
 */
inline static void tinymt64_next_state(tinymt64_t * random) {
    uint64_t x;

    random->status[0] &= TINYMT64_MASK;
    x = random->status[0] ^ random->status[1];
    x ^= x << TINYMT64_SH0;
    x ^= x >> 32;
    x ^= x << 32;
    x ^= x << TINYMT64_SH1;
    random->status[0] = random->status[1];
    random->status[1] = x;
    random->status[0] ^= -((int64_t)(x & 1)) & random->mat1;
    random->status[1] ^= -((int64_t)(x & 1)) & (((uint64_t)random->mat2) << 32);
}

/**
 * This function outputs 64-bit unsigned integer from internal state.
 * Users should not call this function directly.
 * @param random tinymt internal status
 * @return 64-bit unsigned pseudorandom number
 */
inline static uint64_t tinymt64_temper(tinymt64_t * random) {
    uint64_t x;
#if defined(LINEARITY_CHECK)
    x = random->status[0] ^ random->status[1];
#else
    x = random->status[0] + random->status[1];
#endif
    x ^= random->status[0] >> TINYMT64_SH8;
    x ^= -((int64_t)(x & 1)) & random->tmat;
    return x;
}


/**
 * This function certificate the period of 2^127-1.
 * @param random tinymt state vector.
 */
static void period_certification(tinymt64_t * random) {
    if ((random->status[0] & TINYMT64_MASK) == 0 &&
	random->status[1] == 0) {
	random->status[0] = 'T';
	random->status[1] = 'M';
    }
}

/**
 * This function initializes the internal state array with a 64-bit
 * unsigned integer seed.
 * @param random tinymt state vector.
 * @param seed a 64-bit unsigned integer used as a seed.
 */
void tinymt64_init(tinymt64_t * random, uint64_t seed) {
    random->status[0] = seed ^ ((uint64_t)random->mat1 << 32);
    random->status[1] = random->mat2 ^ random->tmat;
    for (int i = 1; i < MIN_LOOP; i++) {
	random->status[i & 1] ^= i + UINT64_C(6364136223846793005)
	    * (random->status[(i - 1) & 1]
	       ^ (random->status[(i - 1) & 1] >> 62));
    }
    period_certification(random);
}

inline static uint64_t tinymt64_generate_uint64(tinymt64_t * random) {
    tinymt64_next_state(random);
    return tinymt64_temper(random);
}

// End TinyMT



// Older Linux kernels don't have CLOCK_MONOTONIC_RAW, so use CLOCK_MONOTONIC instead.
#ifndef CLOCK_MONOTONIC_RAW
#define CLOCK_MONOTONIC_RAW CLOCK_MONOTONIC
#endif

void help()
{
	//       0         1         2         3         4         5         6         7
	//       01234567890123456789012345678901234567890123456789012345678901234567890123456789
	printf ("Usage: hdubench [options] <file>\n");
	printf ("    file is usually a block device, because it doesn't make much sense to\n");
	printf ("    microbenchmark a file on a filesystem.\n");
	printf ("\n");
	printf ("Options: These change settings and define which tests to run\n");
	printf ("    --rpm-measure-time <time>     Default is 3.0. Sets how long to measure the\n");
	printf ("                                  drive RPM in seconds.\n");
	printf ("\n");
	printf ("    --rpm-no-alternate            When measuring RPM, only read sector 0, don't\n");
	printf ("                                  alternate between sectors 0 and 1.\n");
	printf ("\n");
	printf ("    --force-rev-period <time>     Skip RPM measurement and assume the disk\n");
	printf ("                                  revolution time is 'time' microseconds.\n");
	printf ("\n");
	printf ("    --force-sector-size <size>    Assume a different sector size in bytes than\n");
	printf ("                                  what's detected using ioctl(BLKSSZGET)\n");
	printf ("\n");
	printf ("    --skippy <L>                  Run the skippy algorithm with a sequence length\n");
	printf ("                                  of L reads. The skippy algorithm reads sectors\n");
	printf ("                                  0, 1, 3, 6, 10, 15, ... increasing the space\n");
	printf ("                                  between reads linearly.\n");
	printf ("\n");
	printf ("    --seq start,step,end\n");
	printf ("                                  Read an evenly-spaced sequence of sectors.\n");
	printf ("                                  between 'start' and 'end', 'step' sectors apart.\n");
	printf ("\n");
	printf ("    --angpos from,start,step,end,max_error\n");
	printf ("                                  Measure the angular position of each sector\n");
	printf ("                                  relative to 'from'. Measure every 'step' sectors\n");
	printf ("                                  between 'start' and 'end'. Average enough samples\n");
	printf ("                                  so that the standard error of the mean is below\n");
	printf ("                                  'max_error' microseconds. The number of samples\n");
	printf ("                                  is proportional to (1/max_error^2). All 5 numbers\n");
	printf ("                                  must be comma-separated and be a single argument\n");
	printf ("                                  to the --angpos option.\n");
	printf ("\n");
	printf ("    --access from,start,step,end,max_error\n");
	printf ("                                  Measure the access time between sector 'from' and \n");
	printf ("                                  other sectors. This algorithm is nearly the same as\n");
	printf ("                                  angpos, but reports absolute delay instead of\n");
	printf ("                                  just the angular position (mod revolution time)\n");
	printf ("                                  This algorithm is slower due to the inability to\n");
	printf ("                                  measure multiple sectors per disk revolution.\n");
	printf ("                                  Measure every 'step' sectors between 'start'\n");
	printf ("                                  and 'end'. Average enough samples so that the\n");
	printf ("                                  standard error of the mean is below 'max_error'\n");
	printf ("                                  microseconds. The number of samples is proportional\n");
	printf ("                                  to (1/max_error^2). All 5 numbers must be comma-\n");
	printf ("                                  separated and be a single argument to the\n");
	printf ("                                  --access option.\n");
	printf ("\n");
	printf ("    --angpos-list from,max_error\n");
	printf ("                                  Measure the angular position of each sector\n");
	printf ("                                  relative to 'from'. This is the same algorithm as\n");
	printf ("                                  --angpos, but reads a list of sectors to test\n");
	printf ("                                  from stdin, one per line. Each line must start\n");
	printf ("                                  with a number, and the rest of each line is\n");
	printf ("                                  ignored. The output of --track-bounds can\n");
	printf ("                                  be used as the input list to measure the track\n");
	printf ("                                  skew of each track. If 'from' is negative,\n");
	printf ("                                  then the from sector is 'from' sectors before\n");
	printf ("                                  the sector being tested (from=-1 tests track skew).\n");
	printf ("                                  If 'from' is non-negative, it is treated as a constant\n");
	printf ("                                  sector number. Both numbers must be comma-separated\n");
	printf ("                                  and be a single argument to the --angpos-list option.\n");
	printf ("\n");
	printf ("    --access-list from,max_error\n");
	printf ("                                  Measure the access time between sector 'from' and \n");
	printf ("                                  other sectors. This is the same algorithm as\n");
	printf ("                                  --access, but reads a list of sectors to test\n");
	printf ("                                  from stdin, one per line. Each line must start\n");
	printf ("                                  with a number, and the rest of each line is\n");
	printf ("                                  ignored. The output of --track-bounds can\n");
	printf ("                                  be used as the input list to measure the track\n");
	printf ("                                  skew of each track. If 'from' is negative,\n");
	printf ("                                  then the from sector is 'from' sectors before\n");
	printf ("                                  the sector being tested (from=-1 tests track skew).\n");
	printf ("                                  If 'from' is non-negative, it is treated as a constant\n");
	printf ("                                  sector number. Both numbers must be comma-separated\n");
	printf ("                                  and be a single argument to the --access-list option.\n");
	printf ("\n");
	printf ("    --seek-track from,start,step,end\n");
	printf ("                                  Measure the minimum seek time from sector 'from' to\n");
	printf ("                                  a track near another sector. This algorithm tries to\n");
	printf ("                                  find the track-to-track seek time by removing the\n");
	printf ("                                  rotational latency by trying to find the minimum\n");
	printf ("                                  time to access any sector in a region within two\n");
	printf ("                                  tracks of the target sector. 'step' should usually\n");
	printf ("                                  be at least roughly the track size. All 4 numbers\n");
	printf ("                                  must be comma-separated and be a single argument to\n");
	printf ("                                  the --seek-track option.\n");
	printf ("\n");
	printf ("    --track-bounds start,end\n");
	printf ("                                  Find all of the track boundaries and sizes from sector\n");
	printf ("                                  'start' to 'end'. Both numbers must be comma-separated\n");
	printf ("                                  and be a single argument to the --track-bounds option.\n");
	printf ("\n");
	printf ("    --track-bounds-fast start,end\n");
	printf ("                                  A faster version of the track-bounds algorithm. Will skip\n");
	printf ("                                  multiple tracks at a time if the track size hasn't changed.\n");
	printf ("\n");
	printf ("    --random-access start,end,num\n");
	printf ("                                  Do 'num' random access between sectors 'start' to 'end'\n");
	printf ("                                  All numbers must be comma-separated and be a single\n");
	printf ("                                  argument to the --random-access option.\n");
	printf ("\n");

}

class SampleStat {
	double sum, sumsq;
	unsigned int count;

	public:
		SampleStat() {
			reset();
		}
		void reset() {
			sum = 0;
			sumsq = 0;
			count = 0;
		}
		void sample(double value) {
			sum += value;
			sumsq += (value*value);
			count++;
		}
		unsigned int n() {
			return count;
		}
		double mean() {
			return sum/count;
		}
		double stdevp() {
			return sqrt(sumsq/count - (sum/count)*(sum/count));
		}
		double stdev() {
			return sqrt(sumsq/(count-1) - (sum*sum)/(count*(count-1)));
		}
		double stderr() const {
			return sqrt((sumsq - sum*sum/count) / (count*(count-1)));
		}
};

inline unsigned long long get_time_ns()
{
	struct timespec t;
	clock_gettime(CLOCK_MONOTONIC_RAW, &t);
	return t.tv_sec*1000000000ULL + t.tv_nsec;
}
inline unsigned long long get_time_resolution_ns()
{
	struct timespec t;
	clock_getres(CLOCK_MONOTONIC_RAW, &t);
	return t.tv_sec*1000000000ULL + t.tv_nsec;
}
unsigned long long time_ns_pread(int fd, void *buf, size_t nbyte, off_t offset)
{
	struct timespec st, et;
	clock_gettime(CLOCK_MONOTONIC_RAW, &st);

	if (-1 == pread(fd, buf, nbyte, offset)) {
		perror("pread failed");
		return -1;
	}
	clock_gettime(CLOCK_MONOTONIC_RAW, &et);
	return (et.tv_sec - st.tv_sec)*1000000000ULL + (et.tv_nsec - st.tv_nsec);
}
unsigned long long time_ns_pread_abs(int fd, void *buf, size_t nbyte, off_t offset)
{
	struct timespec et;

	if (-1 == pread(fd, buf, nbyte, offset)) {
		perror("pread failed");
		return -1;
	}
	clock_gettime(CLOCK_MONOTONIC_RAW, &et);
	return et.tv_sec*1000000000ULL + et.tv_nsec;
}



unsigned long long measure_rev_period(int fd, void *buf, unsigned int sector_size, unsigned long long measure_time_ns, bool alternate)
{
	printf ("Measuring RPM for %.1f seconds...\n", measure_time_ns/1e9);

	if (-1 == pread(fd, buf, sector_size, 0)) {
		perror("pread failed");
		return -1;
	}

	unsigned int alt_pos = alternate ? sector_size : 0;
	unsigned long long st, et, mt, mintime=0, maxtime=0;
	unsigned int its = 0;
	mt = et = st = get_time_ns();
	for (its=0; st + measure_time_ns >= et; its++) {
		// Maxtor 7405AV: Can't turn off look-ahead, but looks like simply requesting a different sector will not be served from the cache.
		et = time_ns_pread_abs(fd, buf, sector_size, (its&1) ? alt_pos : 0);
		if (its==0 || et-mt > maxtime) maxtime = et-mt;
		if (its==0 || et-mt < mintime) mintime = et-mt;
		mt = et;
	}
	double rpm = its *60.0 / ((et-st)/1e9);
	printf ("    RPM: %.3f (%.3f us per revolution, %u revolutions in %f seconds, max %.3f, min %.3f)\n", rpm, (et-st)/its/1e3, its, (et-st)/1e9, maxtime*1e-3, mintime*1e-3 );
	if (rpm > 20000) {
		// Seagate 15K.7: RCD and DRA are sufficient to disable read cache and prefetch
		// WD S25: DPTL=0 is needed to stop prefetch. It doesn't seem to obey DRA.
		printf ("RPM seems too high for a hard drive. Make sure disk read caching is disabled (Use hdparm -A0 -a0 /dev/<disk> (ATA/SATA) or sdparm -s RCD,DRA,DPTL=0 /dev/<disk>) (SCSI/SAS)\n");
	}
	else {
		if (maxtime > (et-st)*1.5/its)
			printf ("Maximum revolution time is far above the average. RPM measurement is inaccurate because some read requests took more than one revolution to complete. This can happen if the disk became busy during the measurement.\n");
		if (mintime < (et-st)*0.8/its)
			printf ("Minimum revolution time is far belom the average. RPM measurement may be inaccurate.\n");
	}
	return (et-st)/its + 0.5;
}

void skippy(int fd, void *buf, unsigned int size, unsigned long long L)
{
	printf ("Skippy, L = %llu\n", L);
	unsigned long long pos = 0;

	unsigned int *times = (unsigned int*)malloc(L*sizeof(unsigned int));

	// Init
	for (unsigned int i = 0; i < 10; i++)
		time_ns_pread(fd,buf,size,0);

	// Measure
	for (unsigned int i = 0; i < L; i++) {
		times[i] = time_ns_pread(fd, buf, size, pos);
		pos += i*size;
	}

	for (unsigned int i=0;i<L;i++)
		printf ("%u\t%u\n", i, times[i]/1000);

	free (times);
}

void seq(int fd, void *buf, unsigned int size, unsigned long long start, unsigned int step, unsigned long long end)
{
	// Print out results in chunks of lines to reduce the number of prints that could interfere with timing the disk reads.
	const unsigned int CHUNK = 32;

	unsigned int *times = (unsigned int*)malloc(CHUNK*sizeof(unsigned int));

	for (unsigned long long i=0; start+i*step < end; i+= CHUNK)
	{
		// Init
		time_ns_pread(fd,buf,size,i ? (start+(i-1)*step)*size : 0);
		unsigned long long pos = (start+i*step) * size;
		// Measure
		for (unsigned int j = 0; j<CHUNK && pos < end*size; j++) {
			times[j] = time_ns_pread(fd, buf, size, pos);
			pos += step*size;
		}

		for (unsigned int j=0;j<CHUNK && start+(j+i)*step < end;j++)
			printf ("%llu\t%u\n", start+(j+i)*step, times[j]/1000);
		fflush(stdout);
	}

	free (times);
}



typedef struct {
	SampleStat s;
	SampleStat s_nodelta;
	int first;
	int min;
} angpos_data_t;


//	readsize should normally be set to be equal to size (the sector size). Set readsize bigger if you want to
//	read more than one sector at once. Ensure buf is big enough for the entire read.
void angpos(int fd, void *buf, unsigned int size, unsigned int readsize, unsigned long long jump_from, unsigned long long start, unsigned int step, unsigned long long end, double max_error, bool measure_absolute_time, double revtime, bool suppress_header=false)
{
	if (!suppress_header) {
		printf ("Seek time, %s. Jump from %llu, to (%llu - %llu, step %u) with max error %.2f us\n",
			measure_absolute_time ? "absolute" : "relative to angular position",
			jump_from, start, end, step, max_error);
	}

	const unsigned int BUFSZ = 13;
	angpos_data_t data[1<<BUFSZ];		// A circular queue of results for future sectors. Opportunistically sample several sectors over two revolutions to improve speed.
	memset(data, 0, sizeof(data));
	unsigned int pdata=0;
	jump_from *= size;			// Translate sector number to byte position.
	max_error *= 1000.0;		// We use nanoseconds, user uses microseconds

	unsigned int secondary_jump_size = step ? 50/step+16 : 66;		// In units of step*sector_size
	unsigned long long pass_time_limit = revtime * 1.8;	// Time limit for a pass. We want to be able to return to
													// jump_from for the next pass without waiting for another revolution.
	unsigned long long pass_prev_starttime = 0;

	for (unsigned long long i = 0; start+i*step < end; i++)
	{
		double cur_error = 0;
		bool prev_pass_time_limit_exempt = true;
		int max_retries = 2345;		// Try to get error below max_error, but give up after max_retries.
		do
		{
			unsigned long long pos = (start+i*step);
			pread(fd, buf, readsize, jump_from);		// Do two reads of the reference sector? Ensure we start timing with heads stationary. Unsure if this amount of paranoia is really necessary.
			unsigned long long starttime = time_ns_pread_abs(fd, buf, readsize, jump_from);
			if (!prev_pass_time_limit_exempt && starttime - pass_prev_starttime >= 2.5*revtime) {
				pass_time_limit *= 0.95;
				if (pass_time_limit < revtime)
					pass_time_limit = revtime;
				//printf ("pass: Missed start of next pass: new pass_time_limit = %u. Time was %llu\n", pass_time_limit, starttime - pass_prev_starttime);
			}
			else {
				pass_time_limit += 10000;
				if (pass_time_limit > 2*revtime)
					pass_time_limit = 2*revtime;
			}
			pass_prev_starttime = starttime;
			unsigned long long endtime = starttime;
			unsigned int data_index = pdata;
			unsigned int sectors_this_pass = 0;
			prev_pass_time_limit_exempt = false;
			do {
				unsigned long long t = time_ns_pread_abs(fd, buf, readsize, pos*size);
				if ((endtime != starttime) && (t - endtime > 0.98*revtime)) {
					secondary_jump_size *= 1.05;
					prev_pass_time_limit_exempt = true;
				}
				endtime = t;
				angpos_data_t &d = data[data_index];

				long long delta_time = endtime - starttime;

				d.s_nodelta.sample(delta_time);
				if (d.s.n() == 0) {
					d.first = measure_absolute_time ? delta_time : fmod(delta_time, revtime);
					d.s.sample(0);
					d.min = delta_time;
				}
				else {
					double timemod = remainder(delta_time - (long long)d.first, revtime);
					d.s.sample(timemod);
					if (delta_time < d.min)
						d.min = delta_time;
				}
				const unsigned int scaled_secondary_jump_size = secondary_jump_size >> 3;
				pos += scaled_secondary_jump_size*step;
				data_index = (data_index + scaled_secondary_jump_size) & ((1<<BUFSZ)-1);
				sectors_this_pass += scaled_secondary_jump_size;

				if (sectors_this_pass >= (1<<BUFSZ)) 	// Don't overflow the data array
					break;
			} while ( !measure_absolute_time && (endtime - starttime) < pass_time_limit && pos < end);
			if (secondary_jump_size > 32)
				secondary_jump_size--;
			cur_error = (data[pdata].s.n() > 2) ? data[pdata].s.stderr() : max_error;
		} while (--max_retries && (max_error > 0 && cur_error >= max_error));

		angpos_data_t &d = data[pdata];
		if (measure_absolute_time) {
			//double adj = round(((int)d.min - (int)d.first) / revtime) * revtime;
			double adj = round(((int)d.s_nodelta.mean() - d.first) / revtime) * revtime;
			d.first += adj;
		}
		printf ("%llu\t%.0f\t%.1f\t%u\t%u\t%.0f\t%.0f\n", start+(i*step), (d.first + d.s.mean())/1000.0, d.s.stderr()/1000.0, d.s.n(), data[(pdata+1)&((1<<BUFSZ)-1)].s.n(), d.min*1.e-3, d.s_nodelta.mean()*1e-3);

		memset(&data[pdata], 0, sizeof(angpos_data_t));
		pdata = (pdata + 1) & ((1<<BUFSZ)-1);
		fflush(stdout);
	}
}


void seek_profile(int fd, void *buf, unsigned int size, unsigned long long jump_from, unsigned long long from_sector, unsigned int step, unsigned long long end)
{
	const int blind_jump_size = 64;	// 64 sectors is a good value for modern disks, but this must be smaller than a track.
										// The algoritm makes blind backward jumps of this size from somewhere near the beginning of a track
										// and expects to land no further than somewhere in the previous track.
										// The algorithm searches a nearby decreasing region of sectors to find a local minimum access time.
										// The blind_jump_size is used to jump to the previous track once the "beginning" of the current track is found. 

	printf ("Seek track profile. Jump from %llu, to (%llu - %llu, step %u).\n",
		jump_from, from_sector, end, step);

	for (unsigned int i = 0; (from_sector+i*step) < end; i++) {
		unsigned int outer_min_time = ~0U;
		unsigned long long outer_min_pos = 0;

		for (unsigned int k=0; k<2; k++) {
			unsigned long long target;
			if (k == 0)
				target = ((unsigned long long)i * step + from_sector) * size;
			else if (outer_min_pos > blind_jump_size*size) {
				target = outer_min_pos - blind_jump_size*size;		// Jump some sectors, hopefully to preceding track.
			}
			else
				break;

			unsigned int search_step_size = blind_jump_size;
			unsigned int min_time = 4293000000U;
			unsigned long long min_pos = target;
			const unsigned int retries = 5;	// ~Number of retry attempts
			unsigned int j = retries;

			while (j) {
				// Init
				time_ns_pread(fd, buf, size, jump_from * size);

				// Measure
				unsigned int t = time_ns_pread(fd, buf, size, target);

				if (t < min_time + 1000000) {
					min_pos = target;
					if (t < min_time) {
						min_time = t;
						j = retries;
					}

					if ((from_sector+search_step_size)*size > target) {
						target = from_sector*size;
						if (search_step_size > 1)
							search_step_size >>= 1;
						else
							j--;
					}
					else
						target -= search_step_size * size;
				}
				else {
					if (search_step_size > 1) {
						search_step_size >>= 1;
					}
					else
						j--;
					if (search_step_size*size > min_pos)
						target = min_pos - search_step_size*size;
					else
						target = min_pos;
				}
			}

			if (min_time < outer_min_time) {
				outer_min_time = min_time;
				outer_min_pos = min_pos;
			}
		}

		printf ("%llu\t%u\n", outer_min_pos/size, outer_min_time/1000);
		fflush(stdout);
	}
}

unsigned long long find_next_track_boundary (int fd, void *buf, unsigned int size, unsigned long long lb1, unsigned long long ub1, unsigned int min_step, unsigned int min_step_time, double revtime)
{
	unsigned long long ub = ub1, lb = lb1;
	unsigned int retries = 5;

	while (ub - lb > 1 && retries) {
		if (ub - lb > min_step) {
			unsigned long long next_track_partition1 = 0;
			unsigned long long next_track_partition2 = 0;
			unsigned int it_limit = 64;
			do {
				next_track_partition1 = next_track_partition2;
				//printf ("Finding partition: %u %llu-%llu %u %u\n", next_track_partition1, lb, ub, min_step, min_step_time);
				for (unsigned int j = 0; lb+min_step*j < ub+min_step; j++) {
					unsigned long long s = lb+j*min_step;
					unsigned int t = time_ns_pread(fd, buf, size, s*size);

					if (j == 0) continue;
					if (t > 1.7 * revtime) {
						next_track_partition2 = next_track_partition1+1;	// Force fail
						break;		// Try again...
					}
					else {
						if (t > revtime + min_step_time*0.55) t -= revtime;
						if (t > revtime + min_step_time*0.55) {
							min_step *= 1.02;		// If min_step suddenly becomes too small for this region, try slowly increasing it.
							min_step_time *= 1.02;
							next_track_partition2 = next_track_partition1+1;	// Force fail
							break;		// Try again...
						}
						else if (t > (min_step_time + (0.3+0.1) * revtime)) {		// step + skew(0.3rev) + guardband(0.1rev)
							next_track_partition2 = s;
							break;
						}
						else if (t > min_step_time) {
							next_track_partition2 = s;
							break;
						}
					}
				}
			} while (next_track_partition1 != next_track_partition2 && --it_limit);
			if (next_track_partition2 == 0) return -1;		// No boundary found.
			lb = next_track_partition2-min_step;
			ub = next_track_partition2;
		}
		else {
			// ub-lb is never bigger than min_step.
			unsigned long long mp = (ub + lb) >> 1;
			unsigned int tl, tu;
			if (mp >= min_step + lb1) {
				pread(fd, buf, size, (mp-min_step)*size);		// This is always below lb
				tl = time_ns_pread(fd, buf, size, mp*size);
				tu = time_ns_pread(fd, buf, size, (mp+min_step)*size);	// This is always above ub
				if (tl > 1.7 * revtime || tu > 1.7*revtime) {		// Assumption that a correct read does not take more than this amount of time.
					//printf ("Bad: Took too long, trying again.\n");
					retries--;
				}
				else {
					if (retries < 3) {
						// Desperation: if min_step became too small, it could cause tu/tl to exceed one revolution.
						if (tu > revtime) tu -= revtime;
						if (tl > revtime) tl -= revtime;
					}
					if (tu < min_step_time && tl > min_step_time) {
						ub = mp;		// time < min_step_time means no track boundary in that partition
					}
					else if (tl < min_step_time && tu > min_step_time) {
						lb = mp;
					}
					else if (tu < min_step_time && tl < min_step_time) {
						retries = 0;		// Confused: Both were fast, so there is no track boundary at all?
						//printf ("fast: %u - %u: %u %u, %u\n", lb, ub, tl, tu, min_step_time);
						break;
					}
					else
					{
						// Can become confused if the skew is too small and a read across the track
						// boundary causes a full revolution instead of a decrease in latency.
						if (retries < 3)
							printf ("Confused, trying again %u: (%llu - %llu) %u %u\n", retries, lb, ub, tl, tu);
						lb = lb1;
						ub = ub1;
						retries--;
					}
				}
			}
			else {
				pread(fd, buf, size, lb*size);
				tl = time_ns_pread(fd, buf, size, mp*size);
				tu = time_ns_pread(fd, buf, size, ub*size);
				unsigned int low_time_limit = ((mp-lb)*min_step_time/min_step + 0.47*revtime) * 1.1;		// max_skew + time to read half partition if on same track + guardband

				// If controller/something is too slow, then even two sequential sector accesses that
				// cross a track boundary (including skew) will take > 1 revolution to complete.
				// The completion time should be exactly one rev + track skew. But set the lower bound
				// to 1.05*rev time so we don't misclassify accesses that take exactly one rev + 1 sector
				// (consecutive sectors on the same track). We're assuming that the variation in a measured revolution
				// is less than 5%, and that the track skew is significantly over 5% of a revolution.
				tu = fmod(tu-0.05*revtime, revtime) + 0.05*revtime;
				tl = fmod(tl-0.05*revtime, revtime) + 0.05*revtime;


				if (tu > 0.9 * revtime && tl < low_time_limit) {
					ub = mp;		// One revolution = track boundary isn't in the region we tested.
				}
				else if (tl > 0.9 * revtime && tu < low_time_limit) {
					lb = mp;
				}
				else if (tu > 0.9 * revtime && tl > 0.9 * revtime) {
					retries = 0;
					break;
				}
				else
				{
					// Can become confused if a track switch takes too long (e.g., zone switch)
					// Eventually, the search range will increase enough so that the alternative
					// algorithm (using min_step-spaced reads) can be used, which is more robust
					// to weird track switch times.
					//if (retries < 3) printf ("Confused, trying again %u: (%llu - %llu) %u %u\n", retries, lb, ub, tl, tu);
					lb = lb1;
					ub = ub1;
					retries--;
				}
			}
		}
	}
	if (retries == 0) {
		return -1ULL;
	}
	return ub;
}


// get_min_step: Try to find the smallest step size (in sectors) that can be serviced in the same revolution of the disk,
// with some margin so it occurs reliably. A step size that is too small requires two revolutions to service.
// Two phases: First, find an estimate of min_step by doing pairs of accesses and seeing whether it takes one revolution or two.
//     If two revolutions, increase the estimate (multiply by 1+learning_rate), otherwise decrease the estimate (multiply by (1 - 0.01*learning_rate)
// Second phase: Given min_step, estimate a min_step_time (the time it takes to do two accesses spaced min_step apart), with two standard deviations
//     of margin added. This is used as a threshold elsewhere for detecting whether a pair of accesses min_step apart had any extra space (track skew?) inserted.
unsigned int get_min_step(int fd, void *buf, unsigned int size, unsigned long long start, double revtime, unsigned int *min_step_time)
{
	double min_step = 1;
	double min_step_time_avg = 0;
	double lr = 0.25;
	unsigned long long pos = start*size;
	unsigned long long st = get_time_ns();
	bool phase = 0;
	for (unsigned int n = 0,lim=0; n < 200 && lr > 0.0000001 && lim < 10000;n++) {
		unsigned long long et = time_ns_pread_abs(fd, buf, size, pos);

		unsigned int t = et - st;
		min_step_time_avg = min_step_time_avg * (1.0-lr) + t*lr;
		if (t > 0.8*revtime) {
			if (phase)
				lr *= 0.8;
			phase = 0;
			min_step = min_step * (1.0+lr);
			n=0;
		}
		else {
			min_step = min_step * (1.0-0.01*lr);
			phase = 1;
		}
		pos += ((unsigned int)min_step)*size;
		//printf ("%f: %.2f %.2f  %llu %.1f\n", lr, min_step, min_step_time_avg, t, 0.9*revtime);
		st = et;
	}

	// Add a bit of safety margin
	min_step_time_avg *= ceil(min_step*1.4)/min_step;
	min_step = ceil(min_step*1.4);

	//printf ("Debug: min_step_time_avg = %f\n", min_step_time_avg);

	// Get a more accurate estimate of min_step_time
	SampleStat min_step_stat;
	for (unsigned int i = 0;i < 400;i++) {
		unsigned long long et = time_ns_pread_abs(fd, buf, size, pos);
		pos += ((unsigned int)min_step)*size;
		if (et-st < min_step_time_avg*1.2) {		// Reject samples that look obviously wrong.
			min_step_stat.sample(et-st);
		}
		st = et;
	}
	if (min_step_time) {
		if (min_step_stat.n() > 200) {
			*min_step_time = (unsigned int)(min_step_stat.mean() + min_step_stat.stdev()*2.0);
		}
		else {
			printf ("Calibrate min_step: Too many bad samples?\n");
			printf ("  Count: %u, min_step: %.0f, min_step_time_avg: %.0f, mean+2*stdev time: %.0f\n", min_step_stat.n(), min_step, min_step_time_avg, (min_step_stat.mean() + min_step_stat.stdev()*2.0));
			*min_step_time = (unsigned int)min_step_time_avg;
		}
	}
	return min_step;
}

void track_bounds(int fd, void *buf, unsigned int size, unsigned long long start, unsigned long long end, double revtime, bool fastmode)
{
	// fastmode = 1 is similar to MIMD Zone-finding algorithm: http://www.esos.hanyang.ac.kr/files/publication/journals/international/a6-gim.pdf
	// -- Needs extra paranoia because track sizes can change frequently and the MIMD algorithm can produce incorrect answers.
	// fastmode = 0 searches for the next track boundary. It searches one track at a time.
	//   This algorithm is a O(lg n) algorithm for finding track boundaries, with a prediction of the expected track size to make n small.
	//   If the prediction is wrong, exponentially expand the search window until the track boundary is found.
	// fastmode = 1 is the MIMD algorithm with some extra checks added for better reliability.
	//   It can search multiple tracks ahead, assuming that if both N*track_size and (N-1)*track_size are both
	//   track boundaries, then it is highly likely that there are exactly N tracks with the same track_size in this region.
	//   If the track_size estimate needs to change, it falls back to the baseline fastmode=0 algorithm to compute the next track size.
	printf ("Track boundaries%s. (%llu - %llu)\n", fastmode? " (fast)": "",start, end);

	unsigned long long pos = 0;
	unsigned int track_size;
	unsigned int zone_size_estimate = 1;

	unsigned int min_step = 10;
	unsigned int min_step_time = 0;

	// Find minimum step that doesn't require a complete revolution, then add a guardband.
	for (int calibrate_retries = 3; calibrate_retries >= 0; calibrate_retries--) {
		min_step = get_min_step (fd, buf, size, pos, revtime, &min_step_time);
		if (min_step_time < revtime / 2) break;
		if (calibrate_retries == 0) {
			printf ("Calibrate min_step failed. Too big: %u, t=%u\n", min_step, min_step_time);
			return;
		}
	}

	printf ("Using min_step %u, t=%u\n", min_step, min_step_time);

	pos = start;
	track_size = 100;		// Initial estimate of track size, doesn't have to be particularly accurate. Should be less than 2x the real track size, or this algorithm might find track pairs instead of tracks.
	unsigned int c=0;
	while(pos < end) {

		if (fastmode) {
			unsigned long long new_pos;
			bool allow_increase = true;
			unsigned int zone_size = 0;
			unsigned int track_step = zone_size_estimate * 0.9;
			bool first_iteration = true;
			if (track_step < 4) track_step = 4;
			while (true) {
				const unsigned long long target = pos + track_step*track_size;
				if (target < end) {
					// Large jumps: Check two adjacent tracks bounds to avoid an incorrect track_size guess from
					// being accepted as correct if track_step * wrong_track_size happens to be a multiple of correct_track_size.
					// This is less likely to happen for small jumps, because it requires a bigger error in the track size guess.
					if (track_step > 4) {
						new_pos = find_next_track_boundary(fd, buf, size, target-track_size-1, target-track_size+1, min_step, min_step_time, revtime);
					}
					if (track_step <= 4 || new_pos == target-track_size)
						new_pos = find_next_track_boundary(fd, buf, size, target-1, target+1, min_step, min_step_time, revtime);
					else
						new_pos = -1ULL;
				}
				else
					new_pos = -1ULL;

				if (new_pos == target) {		// The prediction was exactly correct
					for (unsigned long long p = pos + track_size; p <= new_pos; p+= track_size)
							printf ("%llu\t%u\n", p, track_size);
					zone_size += track_step;
					if (allow_increase)
						track_step <<= 1;
					else {
						allow_increase = true;
					}
					pos = new_pos;
					if (first_iteration) {
						// Performance optimization: If the zone estimate (new_pos == target) is good, there are likely very few tracks left in the zone,
						// so shrink the track_step estimate to save one (likely-failing) query.
						track_step = 4;
					}
				}
				else {
					if (track_step == 1)
						break;		// Fast mode fail, go back to finding one track boundary.
					track_step >>= 2;
					if (track_step < 1) track_step = 1;
					allow_increase = false;
				}
				first_iteration = false;
				// printf ("track_step = %u\n", track_step);
			}
			if (zone_size > zone_size_estimate * 1.04 + 10)
				zone_size_estimate = zone_size_estimate * 1.04 + 10;
			else if (zone_size > zone_size_estimate)
				zone_size_estimate = (3*zone_size_estimate + zone_size) >> 2;
			else if (zone_size > zone_size_estimate * 0.8)
				zone_size_estimate = zone_size;
			else if (zone_size > 0)
				zone_size_estimate *= 0.8;
		}

		// change_retries: Keep retrying up to a few times if the track boundary keeps changing the track size.
		// Measure the track boundary at least twice, and accept it only if we get the same answer twice in a row.
		// Normally, a track size change results in only one retry to confirm the change.
		// The purpose is to reduce the chance that random noise can make us choose an incorrect
		// position for the next track boundary. It would be very expensive to retry and confirm
		// every track boundary, so only confirm those that change the track size.
		// Track sizes don't tend to change frequently, so I'm more skeptical when it changes than when it doesn't.
		int change_retries = 10;
		unsigned int new_track_size = track_size;
		unsigned long long new_pos = -1ULL;
		while (change_retries--) {
			unsigned int tolerance = 1;
			while (new_pos == -1ULL) {
				if (tolerance >= new_track_size)
					break;
				new_pos = find_next_track_boundary(fd, buf, size, pos+new_track_size-tolerance, pos+new_track_size+tolerance, min_step, min_step_time, revtime);
				//if (new_pos != -1ULL)
				//	printf ("  found: tolerance = %u, diff = %d\n", tolerance, new_pos - (pos + new_track_size));
				tolerance *= 4;
			}
			// We're desperate. Try 10 times before giving up.
			for (int i=0;i<10 && new_pos == -1ULL;i++) {
				if (i > 1) {
					// Recalibrate min_step too?
					unsigned int new_min_step_time;
					unsigned int new_min_step = get_min_step (fd, buf, size, pos, revtime, &new_min_step_time);
					if (new_min_step_time < revtime/2) {
						min_step = new_min_step;
						min_step_time = new_min_step_time;
						if (i > 2) {
							printf ("Desperate: Increasing min_step_time by %.0f%%\n", (i-2)*7.0);
							min_step_time *= 1+(i-2)*0.07;
						}
						printf ("Recalibrate: Using min_step %u, t=%u\n", min_step, min_step_time);
					}
					else
						printf ("Recalibrate attempt failed: Continuing to use min_step %u, t=%u. Attempt was %u, t=%u\n", min_step, min_step_time, new_min_step, new_min_step_time);				
				}
				// The next track boundary must be between 1 and 20000 sectors away. Practically unconstrained because no disk has anywhere near 20K sectors per track. Yet.
				new_pos = find_next_track_boundary(fd, buf, size, pos+1, pos+20000, min_step, min_step_time, revtime);
			}

			if (new_pos == -1ULL)		// Boundary not found despite 10 retries. Skip ahead and continue. There will be incorrect answers here.
				break;
			if ((new_pos - pos) == new_track_size)	// Boundary found. Track size didn't change or was the same as the last retry attempt. Done.
				break;
			new_track_size = new_pos - pos;	// Boundary found, track size changed, so record the tentative new track size and try again to confirm.
		}

		if (new_pos == -1ULL) {
			pos += track_size;
			printf ("Confused: Can't find next track. Giving up, skipping ahead to sector %llu, and trying again.\n", pos);
			continue;
		}

		// Limit growth of track size estimate. If we make an error on one track, this error can
		// persist indefinitely if the estimated track size is an integer multiple of the real track size.

		if ((new_pos - pos) > track_size * 1.2 + 10)
			track_size = track_size * 1.2 + 10;
		else if ((new_pos - pos) < track_size * 0.9)
			track_size = track_size * 0.9;
		else
			track_size = new_pos - pos;

		if (new_pos > end)
			new_pos = end;		// find_next_track_boundary can overrun the end of disk, so clamp the result then quit.
		printf ("%llu\t%lld\n", new_pos, (new_pos - pos));
		if (c++ == 32) {
			c = 0;
			fflush(stdout);
		}
		pos = new_pos;
	}

	return;
}


// A simple test that just does random (uniformly distributed) accesses to a region of the disk.
void random_access(int fd, void *buf, unsigned int size, unsigned long long start, unsigned long long end, unsigned long long num_accesses)
{
	// Random accesses between [start, end) sectors. Includes start, excludes end.
	if (end <= start) {
		printf ("Random access: end sector must be greater than start sector\n");
		return;
	}

	printf ("Random access: %llu access, sectors %llu to %llu\n", num_accesses, start, end);
	fflush(stdout);
	unsigned long long pos = 0;
	tinymt64_t rnd;
	tinymt64_init(&rnd, 0);

	unsigned long long min_t = ~0ULL, max_t = 0;
	unsigned long long st = get_time_ns();
	unsigned long long t = st;
	for (unsigned long long i = 0; i < num_accesses; i++) {
		pos = (tinymt64_generate_uint64(&rnd) % (end-start) + start) * size;
		unsigned long long t2 = time_ns_pread_abs(fd, buf, size, pos);
		if (min_t > t2 - t)
			min_t = t2 - t;
		if (max_t < t2 - t)
			max_t = t2 - t;
		t = t2;
	}
	unsigned long long et = get_time_ns();
	printf ("%.1f\t%llu\t%llu\tmin %.1f\tmax %.1f\t%llu reads in %.1f us\n", (et-st)*1e-3/num_accesses, start, end, min_t*1e-3, max_t*1e-3, num_accesses, (et-st)*1e-3);
}

typedef enum {
	OP_SKIPPY,
	OP_ANGPOS,
	OP_ANGPOS_LIST,
	OP_SEEK_TRACK,
	OP_TRACK_BOUNDS,
	OP_RANDOM_ACCESS,
	OP_SEQ
} cmd_op_e;

struct cmd {
	cmd_op_e op;

	union {
		struct {
			unsigned long long L;
		} skippy;
		struct {
			// All positions are sector numbers
			unsigned long long jump_from;
			unsigned long long start;
			unsigned int step;
			unsigned long long end;
			double max_error;
			bool absolute_time;
		} angpos;
		struct {
			long long jump_from;
			double max_error;
			bool absolute_time;
		} angpos_list;
		struct {
			unsigned long long jump_from;
			unsigned long long start;
			unsigned int step;
			unsigned long long end;
		} seek_track;
		struct {
			unsigned long long start;
			unsigned int step;
			unsigned long long end;
		} seq;
		struct {
			unsigned long long start;
			unsigned long long end;
			bool fastmode;
		} track_bounds;
		struct {
			unsigned long long start;
			unsigned long long end;
			unsigned long long num_accesses;
		} random_access;
	} args;
} cmds[1024];		// Don't try to run more than this many tests

#define CLAMPU(N,LIMIT) do { if ((N)>(LIMIT)) (N)=(LIMIT); } while(0)
int main(int argc, char *argv[])
{
	const char *fname = NULL;
	unsigned long long measure_rpm_time = 3000000000;		// 3 second default
	bool rpm_alternate = true;
	unsigned int force_sector_size = 0;
	unsigned int sector_size=0;
	char *sector_buf=NULL;
	double revtime = 0;

	unsigned int cmd_len = 0;

	while (1) {
		int option_index = 0;
		static struct option long_options[] = {
			{"seq",  required_argument, 0, 0},
			{"angpos",  required_argument, 0, 0},
			{"angpos-list",  required_argument, 0, 0},
			{"access", required_argument, 0, 0},
			{"access-list", required_argument, 0, 0},
			{"seek-track", required_argument, 0, 0},
			{"track-bounds", required_argument, 0, 0},
			{"track-bounds-fast", required_argument, 0, 0},
			{"skippy", required_argument, 0, 0},
			{"random-access", required_argument, 0, 0},
			{"rpm-measure-time", required_argument, 0, 0},
			{"rpm-no-alternate", no_argument, 0, 0},
			{"force-rev-period", required_argument, 0, 0},
			{"force-sector-size", required_argument, 0, 0},
			{"help", no_argument, 0, 0},
			{0, 0, 0, 0}
		};

		int c = getopt_long(argc, argv, "", long_options, &option_index);
		if (c == -1) break;

		switch(c) {
			case 0:
				if (!strcmp(long_options[option_index].name, "skippy")) {
					cmds[cmd_len].op = OP_SKIPPY;
					if (!optarg || 1 != sscanf (optarg, "%llu",
						&cmds[cmd_len].args.skippy.L)) {
						printf ("Error in arguments to --skippy \"%s\". Expecting 1 number\n", optarg);
						return -2;
					}
					else
						cmd_len++;
				}
				else if (!strcmp(long_options[option_index].name, "seq")) {
					cmds[cmd_len].op = OP_SEQ;
					if (!optarg || 3 != sscanf (optarg, "%llu, %u, %llu",
						&cmds[cmd_len].args.seq.start,
						&cmds[cmd_len].args.seq.step,
						&cmds[cmd_len].args.seq.end)) {
						printf ("Error in arguments to --seq \"%s\". Expecting comma-separated list of 3 numbers\n", optarg);
						return -2;
					}
					else {
						cmd_len++;
					}
				}
				else if (!strcmp(long_options[option_index].name, "angpos") ||
						!strcmp(long_options[option_index].name, "access")) {
					cmds[cmd_len].op = OP_ANGPOS;
					cmds[cmd_len].args.angpos.absolute_time = long_options[option_index].name[1] == 'c';	// aNgpos vs. aCcess
					if (!optarg || 5 != sscanf (optarg, "%llu, %llu, %u, %llu, %lf",
						&cmds[cmd_len].args.angpos.jump_from,
						&cmds[cmd_len].args.angpos.start,
						&cmds[cmd_len].args.angpos.step,
						&cmds[cmd_len].args.angpos.end,
						&cmds[cmd_len].args.angpos.max_error)) {
						printf ("Error in arguments to --%s \"%s\". Expecting comma-separated list of 5 numbers\n", long_options[option_index].name, optarg);
						return -2;
					}
					else {
						cmd_len++;
					}
				}
				else if (!strcmp(long_options[option_index].name, "angpos-list") ||
						!strcmp(long_options[option_index].name, "access-list")) {
					cmds[cmd_len].op = OP_ANGPOS_LIST;
					cmds[cmd_len].args.angpos_list.absolute_time = long_options[option_index].name[1] == 'c';	// aNgpos vs. aCcess
					if (!optarg || 2 != sscanf (optarg, "%lld, %lf",
						&cmds[cmd_len].args.angpos_list.jump_from,
						&cmds[cmd_len].args.angpos_list.max_error)) {
						printf ("Error in arguments to --%s \"%s\". Expecting comma-separated list of 2 numbers\n", long_options[option_index].name, optarg);
						return -2;
					}
					else {
						cmd_len++;
					}
				}
				else if (!strcmp(long_options[option_index].name, "seek-track")) {
					cmds[cmd_len].op = OP_SEEK_TRACK;
					if (!optarg || 4 != sscanf (optarg, "%llu, %llu, %u, %llu",
						&cmds[cmd_len].args.seek_track.jump_from,
						&cmds[cmd_len].args.seek_track.start,
						&cmds[cmd_len].args.seek_track.step,
						&cmds[cmd_len].args.seek_track.end)) {
						printf ("Error in arguments to --seek-track \"%s\". Expecting comma-separated list of 4 numbers\n", optarg);
						return -2;
					}
					else {
						cmd_len++;
					}
				}
				else if (!strcmp(long_options[option_index].name, "track-bounds") ||
						!strcmp(long_options[option_index].name, "track-bounds-fast")) {
					cmds[cmd_len].op = OP_TRACK_BOUNDS;
					cmds[cmd_len].args.track_bounds.fastmode = (long_options[option_index].name[12] != 0);
					if (!optarg || 2 != sscanf (optarg, "%llu, %llu",
						&cmds[cmd_len].args.track_bounds.start,
						&cmds[cmd_len].args.track_bounds.end)) {
						printf ("Error in arguments to %s \"%s\". Expecting comma-separated list of 2 numbers (start sector, end sector)\n", long_options[option_index].name, optarg);
						return -2;
					}
					else {
						cmd_len++;
					}
				}
				else if (!strcmp(long_options[option_index].name, "random-access")) {
					cmds[cmd_len].op = OP_RANDOM_ACCESS;
					if (!optarg || 3 != sscanf (optarg, "%llu, %llu, %llu",
						&cmds[cmd_len].args.random_access.start,
						&cmds[cmd_len].args.random_access.end,
						&cmds[cmd_len].args.random_access.num_accesses)) {
						printf ("Error in arguments to %s \"%s\". Expecting comma-separated list of 3 numbers (start sector, end sector, number of accesses)\n", long_options[option_index].name, optarg);
						return -2;
					}
					else {
						cmd_len++;
					}
				}
				else if (!strcmp(long_options[option_index].name, "rpm-measure-time")) {
					double s;
					if (!optarg || 1 != sscanf (optarg, "%lf", &s)) {
						printf ("Error in arguments to --rpm-measure-time \"%s\". Expecting one floating-point number\n", optarg);
						return -2;
					}
					measure_rpm_time = (s * 1e9) + 0.5;
				}
				else if (!strcmp(long_options[option_index].name, "rpm-no-alternate")) {
					rpm_alternate = false;
				}
				else if (!strcmp(long_options[option_index].name, "force-rev-period")) {
					double s;
					if (!optarg || 1 != sscanf (optarg, "%lf", &s)) {
						printf ("Error in arguments to --force-rev-period \"%s\". Expecting one floating-point number\n", optarg);
						return -2;
					}
					revtime = (s * 1e3);		// We use nanoseconds internally
				}
				else if (!strcmp(long_options[option_index].name, "force-sector-size")) {
					if (!optarg || 1 != sscanf (optarg, "%u", &force_sector_size)) {
						printf ("Error in arguments to --force-sector-size \"%s\". Expecting one unsigned integer\n", optarg);
						return -2;
					}
				}
				else if (!strcmp(long_options[option_index].name, "help")) {
					help();
					return 0;
				}
				else {
					printf ("Bad option %s. Should not happen.\n", long_options[option_index].name);
					return -1;
				}
				break;
			default:
				return -1;
				break;
		}

		if (cmd_len >= (sizeof(cmds)/sizeof(cmds[0]))) {		// Prevent overflow.
			printf ("Too many commands (%d)\n", cmd_len);
			return -1;
		}
	}
	while (optind < argc) {
		fname = argv[optind++];
	}
	if (!fname) {
		printf ("Filename (block device to measure) is required\n");
		return -1;
	}

	printf ("Opening %s\n", fname);

	int fd = open(fname, O_DIRECT | O_RDONLY);
	if (fd == -1)	{
		perror("Failed to open");
		return -1;
	}
	posix_fadvise(fd, 0, 0, POSIX_FADV_RANDOM);

	unsigned long long sz = -1;
	sz = lseek(fd, 0, SEEK_END);
	ioctl(fd, BLKSSZGET, &sector_size);		// This is usually wrong for 4KB sector drives that emulate 512B. Manually force 4K when testing one of those drives.

	printf ("Sector size is %u\n", sector_size);

	if (force_sector_size) {
		printf ("Forcing sector size to be %u bytes\n", force_sector_size);
		sector_size = force_sector_size;
	}
	if (sector_size == 0) {
		printf ("Sector size can't be zero. Use --force-sector-size to choose a sector size if necessary.\n");
		return -1;
	}


	unsigned long long sz_sectors = sz / sector_size;
	printf ("File size is %llu bytes, %llu sectors\n", sz, sz_sectors);
	printf ("Timer claims to have a resolution of %llu ns\n", get_time_resolution_ns());


	// Start by measuring RPM, because it's used by some of the other tests
	sector_buf = (char*)mmap(NULL, sector_size*16, PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS, 0, 0);
	if (revtime == 0) {
		revtime = measure_rev_period(fd, sector_buf, sector_size, measure_rpm_time, rpm_alternate);
	}
	else {
		printf ("Assuming RPM = %.3f (%.3f us per revolution)\n", 60e9/revtime, revtime/1e3);
	}


	printf ("Running %u tests...\n", cmd_len);

	for (unsigned int t = 0; t < cmd_len; t++) {
		unsigned long long test_runtime_start = get_time_ns();
		switch (cmds[t].op) {
			case OP_SKIPPY:
				CLAMPU(cmds[t].args.skippy.L, sz_sectors-1);
				skippy (fd, sector_buf, sector_size, cmds[t].args.skippy.L);
				break;
			case OP_SEQ:
				CLAMPU(cmds[t].args.seq.start, sz_sectors-1);
				CLAMPU(cmds[t].args.seq.end, sz_sectors);
				seq(fd, sector_buf, sector_size,
					cmds[t].args.seq.start,
					cmds[t].args.seq.step,
					cmds[t].args.seq.end);
				break;
			case OP_ANGPOS:
				CLAMPU(cmds[t].args.angpos.jump_from, sz_sectors-1);
				CLAMPU(cmds[t].args.angpos.start, sz_sectors-1);
				CLAMPU(cmds[t].args.angpos.end, sz_sectors);

				angpos(fd, sector_buf, sector_size, sector_size,
					cmds[t].args.angpos.jump_from,
					cmds[t].args.angpos.start,
					cmds[t].args.angpos.step,
					cmds[t].args.angpos.end,
					cmds[t].args.angpos.max_error,
					cmds[t].args.angpos.absolute_time,
					revtime);

				break;
			case OP_ANGPOS_LIST:
				if (cmds[t].args.angpos_list.jump_from >= (long long)sz_sectors)
					cmds[t].args.angpos_list.jump_from = (long long)sz_sectors-1;
				printf ("Access time list, %s. Jump from %lld with max error %.2f us\n",
					cmds[t].args.angpos_list.absolute_time ? "absolute" : "angular position",
					cmds[t].args.angpos_list.jump_from, cmds[t].args.angpos_list.max_error);
				{
					char line[1024];
					while (fgets(line, sizeof(line), stdin)) {
						unsigned long long target;
						if (sscanf(line, "%llu", &target) == 1 && target < sz_sectors) {
							unsigned long long from = cmds[t].args.angpos_list.jump_from;
							if (cmds[t].args.angpos_list.jump_from < 0) from += target;
							if (from < sz_sectors) {
								angpos(fd, sector_buf, sector_size, sector_size,
									from,
									target,
									1,
									target+1,
									cmds[t].args.angpos_list.max_error,
									cmds[t].args.angpos_list.absolute_time,
									revtime, true);		// suppress header message
							}
						}
					}
				}
				break;
			case OP_SEEK_TRACK:
				CLAMPU(cmds[t].args.seek_track.jump_from, sz_sectors-1);
				CLAMPU(cmds[t].args.seek_track.start, sz_sectors-1);
				CLAMPU(cmds[t].args.seek_track.end, sz_sectors);
				seek_profile(fd, sector_buf, sector_size,
					cmds[t].args.seek_track.jump_from,
					cmds[t].args.seek_track.start,
					cmds[t].args.seek_track.step,
					cmds[t].args.seek_track.end);
				break;
			case OP_TRACK_BOUNDS:
				CLAMPU(cmds[t].args.track_bounds.start, sz_sectors-1);
				CLAMPU(cmds[t].args.track_bounds.end, sz_sectors);
				track_bounds(fd, sector_buf, sector_size,
					cmds[t].args.track_bounds.start,
					cmds[t].args.track_bounds.end,
					revtime,
					cmds[t].args.track_bounds.fastmode);
				break;
			case OP_RANDOM_ACCESS:
				CLAMPU(cmds[t].args.random_access.end, sz_sectors);
				CLAMPU(cmds[t].args.random_access.start, cmds[t].args.random_access.end);
				random_access(fd, sector_buf, sector_size,
					cmds[t].args.random_access.start,
					cmds[t].args.random_access.end,
					cmds[t].args.random_access.num_accesses);
				break;
			default:
				printf ("Unknown test %u\n", cmds[t].op);
				break;
		}
		printf ("Test runtime: %.3f s\n", (get_time_ns()-test_runtime_start)/1.0e9);
	}

	close (fd);
	return 0;
}
