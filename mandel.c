/*
    Name: Shivam Patel
    ID: 1001707265
*/

#include "bitmap.h"

#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <pthread.h>
#include <sys/time.h>

// the struct contains all the parameters required for the thread
struct parameters{
    int max;
    int begin;
    int end;
    double xmin;
    double xmax;
    double ymin;
    double ymax;
};

int iteration_to_color( int i, int max );
int iterations_at_point( double x, double y, int max );
void * computing_image( void * arg );
struct bitmap *bm;


void show_help()
{
	printf("Use: mandel [options]\n");
	printf("Where options are:\n");
	printf("-m <max>    The maximum number of iterations per point. (default=1000)\n");
	printf("-x <coord>  X coordinate of image center point. (default=0)\n");
	printf("-y <coord>  Y coordinate of image center point. (default=0)\n");
	printf("-s <scale>  Scale of the image in Mandlebrot coordinates. (default=4)\n");
	printf("-W <pixels> Width of the image in pixels. (default=500)\n");
	printf("-H <pixels> Height of the image in pixels. (default=500)\n");
	printf("-o <file>   Set output file. (default=mandel.bmp)\n");
    printf("-n <no_of_threads>  The number of threads\n");
	printf("-h          Show this help text.\n");
	printf("\nSome examples are:\n");
	printf("mandel -x -0.5 -y -0.5 -s 0.2\n");
	printf("mandel -x -.38 -y -.665 -s .05 -m 100\n");
	printf("mandel -x 0.286932 -y 0.014287 -s .0005 -m 1000\n\n");
}

int main( int argc, char *argv[] )
{
    // struct to calculate the time for each mandelbrot
    struct timeval start;
    struct timeval finish;
	char c;

	// These are the default configuration values used
	// if no command line arguments are given.

	const char *outfile = "mandel.bmp";
	double xcenter = 0;
	double ycenter = 0;
	double scale = 4;
	int    image_width = 500;
	int    image_height = 500;
	int    max = 1000;
    // initially the number of threads is 1
    int no_of_threads = 1; 
    int i = 0;

	// For each command line argument given,
	// override the appropriate configuration value.

	while((c = getopt(argc,argv,"x:y:s:W:H:m:o:n:h"))!=-1) {
		switch(c) {
			case 'x':
				xcenter = atof(optarg);
				break;
			case 'y':
				ycenter = atof(optarg);
				break;
			case 's':
				scale = atof(optarg);
				break;
			case 'W':
				image_width = atoi(optarg);
				break;
			case 'H':
				image_height = atoi(optarg);
				break;
			case 'm':
				max = atoi(optarg);
				break;
			case 'o':
				outfile = optarg;
				break;
            case 'n': // the -n in the argument to get the number of threads
				no_of_threads = atoi(optarg);
				break;
			case 'h':
				show_help();
				exit(1);
				break;
		}
	}

    gettimeofday( &start, NULL );

	// Display the configuration of the image.
    // updated the printf to add no of threads in the output
	printf("mandel: x=%lf y=%lf scale=%lf max=%d Threads=%d outfile=%s\n",xcenter,ycenter,scale,max,no_of_threads,outfile);

	// Create a bitmap of the appropriate size.
	bm = bitmap_create(image_width,image_height);
    
	// Fill it with a dark blue, for debugging
	bitmap_reset(bm,MAKE_RGBA(0,0,255,0));

	// Compute the Mandelbrot image
	//computing_image(bm,xcenter-scale,xcenter+scale,ycenter-scale,ycenter+scale,max,no_of_threads);
    
    struct parameters params[no_of_threads];
    pthread_t tid[no_of_threads];
    // if the number of threads is 1 which by default is also 1
    // it would execute first for loop
    // for threads > 1 it would got the else condition and execute the other for loop
    if(no_of_threads == 1)
    {
        for(i = 0; i < no_of_threads; i++)
        {
            params[i].xmin = xcenter - scale;
            params[i].xmax = xcenter + scale;
            params[i].ymin = ycenter - scale;
            params[i].ymax = ycenter + scale;
            params[i].max = max;
            params[i].begin = 0;
            //printf("begin <= 1: %d \n",h1 * i);
            params[i].end = image_height;
            //printf("begin <= 1: %d \n",h1 * i);
            pthread_create( &tid[i], NULL, computing_image, (void*) &params[i] );
        } 
    }
    else
    {
        for(i = 0; i < no_of_threads; i++)
        {
            // assigning values to the thread arguments
            // for x axis xmin is the leftmost side of the image and while xmax is the righmost side of the image
            // similarly for y axis ymin is lowermost and ymax is topmost
            params[i].xmin = xcenter - scale;
            params[i].xmax = xcenter + scale;
            params[i].ymin = ycenter - scale;
            params[i].ymax = ycenter + scale;
            params[i].max = max;  
            params[i].begin = (image_height / no_of_threads) * i;
            //printf("begin > 1: %d \n",(image_height / no_of_threads) * i);  
            params[i].end = (image_height / no_of_threads) * (i + 1);  
            //printf("end > 1 %d \n",(image_height / no_of_threads) * (i + 1));    
            pthread_create( &tid[i], NULL, computing_image, (void*) &params[i] );
        }
    }
    for( i = 0; i < no_of_threads; i++ )
    { 
        //printf("No of threads %d \n",i);
        pthread_join( tid[i], NULL );
    }
	// Save the image in the stated file.
	if(!bitmap_save(bm,outfile)) {
		fprintf(stderr,"mandel: couldn't write to %s: %s\n",outfile,strerror(errno));
		return 1;
	}

    // get the time it takes for each creation of mandelbrot
    gettimeofday( &finish, NULL );

    int time_duration = ( ( finish.tv_sec - start.tv_sec ) * 1000000 + ( finish.tv_usec - start.tv_usec ) );

    printf("Duration: %d\n", time_duration );

	return 0;
}

/*
Compute an entire Mandelbrot image, writing each point to the given bitmap.
Scale the image to the range (xmin-xmax,ymin-ymax), limiting iterations to "max"
*/

void * computing_image( void * arg )
{
  int i,j;

  struct parameters * params = (struct parameters *) arg;

  int width = bitmap_width(bm);
  int height = bitmap_height(bm);
  int begin = params -> begin;
  int end = params -> end;
  int max = params -> max;
  double xmin = params -> xmin;
  double xmax = params -> xmax;
  double ymin = params -> ymin;
  double ymax = params -> ymax;

  for(j=begin;j<end;j++) {

		for(i=0;i<width;i++) {

			// Determine the point in x,y space for that pixel.
			double x = xmin + i*(xmax-xmin)/width;
			double y = ymin + j*(ymax-ymin)/height;
            //printf("i %d \n",i);
            //printf("j %d \n",j);
			// Compute the iterations at that point.
			int iters = iterations_at_point(x,y,max);
            //printf("iters %d \n",iters);
			// Set the pixel in the bitmap.
			bitmap_set(bm,i,j,iters);
		}
	}

  return NULL;
}


/*
Return the number of iterations at point x, y
in the Mandelbrot space, up to a maximum of max.
*/

int iterations_at_point( double x, double y, int max )
{
	double x0 = x;
	double y0 = y;

	int iter = 0;

	while( (x*x + y*y <= 4) && iter < max ) {

		double xt = x*x - y*y + x0;
		double yt = 2*x*y + y0;

		x = xt;
		y = yt;

		iter++;
	}

	return iteration_to_color(iter,max);
}

/*
Convert a iteration number to an RGBA color.
Here, we just scale to gray with a maximum of imax.
Modify this function to make more interesting colors.
*/

int iteration_to_color( int i, int max )
{
	int gray = 255*i/max;
	return MAKE_RGBA(gray,gray,gray,0);
}
