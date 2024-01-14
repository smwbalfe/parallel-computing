#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

struct Context {
    uint array_size;
    int n_processors;
    int rank;
    double precision;
    int *block_size;
    int *displacements;
    double *local_buffer;
    double *input_buffer;
    int complete;
} typedef context_t;


void print_array(void* input_buffer, uint array_size) {
    double *in = (double *) input_buffer;
    for (uint y = 0; y < array_size; ++y) {
        for (int x = 0; x < array_size; ++x) {
            printf("%f ", in[array_size * y + x]);
        }
        printf("\n");
    }
}

double calculate_average(const double *input_buffer, const int y, const int x, const uint array_size){

    double above = input_buffer[array_size * (y-1) + x];
    double left = input_buffer[array_size * y + (x-1)];
    double below = input_buffer[array_size * (y+1) + x];
    double right = input_buffer[array_size * y + (x+1)];

    return (above + left + below + right) / 4;
}

void array_relaxation(context_t *context){

    int array_offset = context->displacements[context->rank];

    // go through each value in our local buffer
    for (int i = 0; i < context->block_size[context->rank] ; ++i){

        // calculate what position we are at in the initial input buffer.
        int y = (array_offset + i) / (int) context->array_size;
        int x = (array_offset + i) % (int) context->array_size;

        // check for borders
        if (y == 0 || x == 0 || y == context->array_size - 1 || x == context->array_size - 1){
            continue;
        }

        double old_value = context->local_buffer[i];
        double new_value = calculate_average(context->input_buffer,y,x,context->array_size);
        context->local_buffer[i] = new_value;

        // set complete to 0 if the change is greater than the desired precision ,signalling that all threads must complete
        // another cycle.
        if (fabs(new_value - old_value) > context->precision){
            context->complete = 0;
        }
    }
}

int main(int argc, char **argv) {

    int rc = MPI_Init(&argc, &argv);

    if (rc != MPI_SUCCESS) {
        printf("Error starting MPI test program\n");
        MPI_Abort(MPI_COMM_WORLD, rc);
    }

    // create our context which just allows us to easily pass variables to other functions, no more use than that
    context_t *context = malloc(sizeof(context_t));
    context->array_size = (int) strtol(argv[1], 0,10);
    context->precision = strtof(argv[2], 0);

    // store the rank of the processor alongside the total count
    MPI_Comm_rank(MPI_COMM_WORLD, &context->rank);
    MPI_Comm_size(MPI_COMM_WORLD, &context->n_processors);

    // allocate the initial arrays to store rank block size, rank displacement in buffer and initial input buffer.
    context->block_size = malloc((ssize_t) sizeof(double) * context->n_processors);
    context->displacements = malloc((ssize_t) sizeof(double) * context->n_processors);
    context->input_buffer = malloc((ssize_t) sizeof(double) * ((ssize_t) pow(context->array_size,2)));

    // rank 0 calculates displacements and block sizes
    if (context->rank == 0) {

        uint remainder = (uint) (pow(context->array_size,2)) % context->n_processors;

        int sum = 0;

        // go through each processor
        for (uint i = 0; i < context->n_processors; ++i) {

            // block size will be the initial division
            context->block_size[i] = (int) (pow(context->array_size,2)) / context->n_processors;

            // if we have remainder values left, increment by 1
            if (remainder > 0) {
                context->block_size[i]++;
                remainder--;
            }

            // displacements is just the running sum we are at in terms of workload division
            context->displacements[i] = sum;

            // sum is just the block size addition for said thread
            sum += context->block_size[i];
        }
    }

    // broadcast block size and displacements.
    MPI_Bcast(context->block_size, context->n_processors, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(context->displacements, context->n_processors, MPI_INT, 0, MPI_COMM_WORLD);

    context->local_buffer = malloc((ssize_t) sizeof(double) * context->block_size[context->rank]);

    // make one processor allocate the array
    if (context->rank == 0) {
        for (int y = 0; y < context->array_size; ++y) {
            if (y == 0) {
                for (int x = 0; x < context->array_size; ++x) {
                    context->input_buffer[context->array_size * y + x] = 1;
                }
            } else {
                for (int x = 0; x < context->array_size; ++x) {
                    if (x == 0) {
                        context->input_buffer[context->array_size * y + x] = 1;
                    } else {
                        context->input_buffer[context->array_size * y + x] = 2;
                    }
                }
            }
        }
        context->complete = 1;
    }

    MPI_Bcast(&context->complete, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(context->input_buffer,  (int) pow(context->array_size,2), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    int result = 0;
    double start_time = MPI_Wtime();

    // loop forever until break condition met
    while(1) {

        // scatter the input buffer according to the displacement on where to start reading from and the block size of how many elements to read
        // from that position
        MPI_Scatterv(context->input_buffer, (int*) context->block_size, (int*) context->displacements, MPI_DOUBLE,
                     context->local_buffer, (int) context->block_size[context->rank], MPI_DOUBLE,
                     0, MPI_COMM_WORLD);

        // perform an array relaxation pass
        array_relaxation(context);

        // same as scatter simply read the displacements to and block size to figure out the size of the buffer to read into the input buffer.
        MPI_Gatherv(context->local_buffer, context->block_size[context->rank],
                    MPI_DOUBLE, context->input_buffer,
                    (int*) context->block_size, (int*)context->displacements,
                    MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // after each iteration, all threads must have an updated input buffer copy to use on their next iteration
        MPI_Bcast(context->input_buffer, (int) pow(context->array_size,2), MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // Reduce and broadcast, we reduce the complete flag of each processor by summing all the values up
        MPI_Allreduce(&context->complete, &result, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        // if the sum of all complete flags is the number of processors , i.e. every processor has a complete flag to set 1, we can break
        // out as we are complete
        if (result == context->n_processors) {
            break;
        } else {
            // not all processors are complete , reset the flag and do another cycle.
            context->complete = 1;
        }
    }

    double end_time = MPI_Wtime();
    double elapsed = end_time - start_time;

    double total_elapsed_time;

    MPI_Reduce(&elapsed, &total_elapsed_time, 1 , MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (context->rank == 0) {
        print_array(context->input_buffer, context->array_size);
        //double average = total_elapsed_time / context->n_processors;
        //printf("time elapsed: %f\n", average);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}
