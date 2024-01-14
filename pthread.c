#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <math.h>
#include <string.h>
#include <time.h>

struct Cell {
    double value;
} typedef Cell;


struct Context {
    int threadCount_;
    int arraySize_;
    double precision_;
    int* operationsPerThread; // array of how many operations each thread performs
    pthread_t* threads;
    pthread_barrier_t barrier;
    int complete; // controls whether the algorithm is complete
    Cell* inputBuffer; // input buffer for threads to read
    Cell* outputBuffer; // output buffer for threads to write too
    int iterCount;
    int mainThreadComplete;
} typedef Context;


void print_array(void* arr, int arrSize) {

    Cell* array = (Cell *) arr;
    for (int y = 0; y < arrSize; ++y) {

        for (int x = 0; x < arrSize; ++x) {
            printf("%f ", array[arrSize * y + x].value);
        }
        printf("\n");
    }
}


double setAverage(Cell* arr, int y, int x, int arrSize){
    double above = arr[arrSize * (y-1) + x].value;
    double left = arr[arrSize * y + (x-1)].value;
    double below = arr[arrSize * (y+1) + x].value;
    double right = arr[arrSize * y + (x+1)].value;


    return (above + left + below + right) / 4;
}


struct threadContext {
    unsigned int start; // index where the thread should start work on
    unsigned int blockSize; // index where thread should end.
    Context* context; // data used between all threads that is read only
    unsigned int id; // id of the thread if required
} typedef threadContext;

/* given an itemNumber of the inner array returns it original position in the border array

  example = 5 x 5 array including border
  3 x 3 is the 5 x 5 inner array so itemNumber 0 = the first value of 3 x 3
  itemNumber 0 would give the index of the original 5 x 5 full array by using a few jump calculations to ignore borders.
  in this case it would 0 would give 6
 * */
int getInnerArrayIndex(int itemNumber, Context* ctx){
    return (ctx->arraySize_ + 1) + itemNumber + (itemNumber / (ctx->arraySize_ - 2) * 2);

}


void work(void* args) {

    threadContext *threadArgs = (threadContext *) args;

    /* infinite while loop from which threads break out of when they have been signalled to complete by the main thread*/
    while(1) {

        /* start from the specified index position and go blocksize places from there*/
        for (int cell = threadArgs->start; cell < (threadArgs->start + threadArgs->blockSize); cell++) {

            /* convert this value into the index it would be if it were using the original border array*/
            int convertIndex = getInnerArrayIndex(cell, threadArgs->context);

            /* figure out what position we are at in the array in terms of rows and columns */
            int y = convertIndex / threadArgs->context->arraySize_;
            int x = convertIndex % threadArgs->context->arraySize_;

            /* each thread now reads its specified value from teh read only input buffer */
            Cell *current = &threadArgs->context->inputBuffer[convertIndex];

            /* hold the old value for precision check*/
            double old = current->value;


            /* calculate the average of its four neighbours*/
            double new = setAverage(threadArgs->context->inputBuffer, y, x, threadArgs->context->arraySize_);

            /* write to shared output buffer the new value*/
            threadArgs->context->outputBuffer[convertIndex].value = new;

            /* if any thread has a precision that is greater than the specified value then we are not complete*/
            if ((fabs(new-old) > threadArgs->context->precision_)) {
                threadArgs->context->complete = 0;
            }
        }
        /* first barrier ensures that we wait for each thread to finish its iteration before continuing*/
        pthread_barrier_wait(&threadArgs->context->barrier);

        /*  second barrier makes the worker threads wait for the main thread to do the sequential operation
         *
         *  > This sequential operation involves swapping the buffers, resetting the complete flag to 1 and
         *    setting a flag of completion if the complete flag stays at 1 after all threads have completed
         * */

        pthread_barrier_wait(&threadArgs->context->barrier);

        /* if main thread has set this flag while sequentially running, then we signal all thread to exit at the same time*/
        if (threadArgs->context->mainThreadComplete){
            break;
        }

    }
}

void controlThread(Context* context) {
    while (1) {
        /* wait on all worker thread to finish their iteration*/
        pthread_barrier_wait(&context->barrier);

        /* if at least one worker changed the complete flag we reset that flag and swap the buffers around*/
        if (!context->complete){

            context->iterCount++;
            context->complete = 1;

            Cell *temp = context->inputBuffer ;
            context->inputBuffer = context->outputBuffer;
            context->outputBuffer = temp;

            /* all the worker threads will be waiting at this barrier so main signals now to move on*/
            pthread_barrier_wait(&context->barrier);
        }
        else{
            /* if complete flag changed, set a flag to signal completion*/
            context->mainThreadComplete = 1;

            /* wait on worker threads in the work function, they will now see this flag set and break
             * at the same time as the main does.
             * */
            pthread_barrier_wait(&context->barrier);

            /* exit */
            break;
        }
    }
}

int main(int argc, char** argv) {

    /* fill context with thread count, array dimensions and precision required for this run*/
    Context context;

    context.arraySize_ = (int) strtol(argv[2], 0, 10);
    context.threadCount_ = (int) strtol(argv[2], 0, 10);
    context.precision_ = strtof(argv[3], 0);
    context.operationsPerThread = malloc(context.threadCount_ * sizeof(int));
    context.threads = malloc(sizeof(pthread_t) * context.threadCount_);
    context.complete = 0;
    context.mainThreadComplete = 0;

    /* divide the workload up between threads*/

    /* barrier used between the main thread and worker threads so ctx.threadCount_ + 1*/
    pthread_barrier_init(&context.barrier, NULL, context.threadCount_ + 1);

    Cell* inputBuffer = malloc(sizeof(struct Cell) * (context.arraySize_ * context.arraySize_));


    /* initial input buffer filled with border values and set value */

    for (int y = 0; y < context.arraySize_; ++y) {
        if (y == 0) {
            for (int x = 0; x < context.arraySize_; ++x) {
                inputBuffer[context.arraySize_ * y + x].value = 1;
            }
        } else {
            for (int x = 0; x < context.arraySize_; ++x) {
                if (x == 0) {
                    inputBuffer[context.arraySize_ * y + x].value = 1;
                } else {
                    inputBuffer[context.arraySize_ * y + x].value = 0;
                }
            }
        }
    }


    Cell *outputBuffer = (Cell *) malloc(sizeof(struct Cell) * (context.arraySize_ * context.arraySize_));


    memcpy(outputBuffer, inputBuffer, sizeof(struct Cell) * context.arraySize_ * context.arraySize_);
    context.inputBuffer = inputBuffer;
    context.outputBuffer = outputBuffer;

    /* divide the inner array excluding the borders by total number of threads*/
    int ops = ((context.arraySize_ - 2) * (context.arraySize_ - 2)) / context.threadCount_;

    /* obtain the leftover values for when it does not perfectly divide */
    int remainder = ((context.arraySize_ - 2) * (context.arraySize_ - 2)) % context.threadCount_;

    /* assign work to each thread*/
    for (int i = 0; i < context.threadCount_; ++i) {
        /* if there are some remainder values left then add an extra operation for the thread to work on
         * and decrement remainder until there are non left*/
        if (remainder != 0) {
            context.operationsPerThread[i] = ops + 1;
            remainder--;
        } else {
            context.operationsPerThread[i] = ops;
        }
    }
    struct timespec start, finish;
    double elapsed;
    clock_t CPUtime;

    clock_gettime(CLOCK_MONOTONIC, &start);
    CPUtime = clock();

    /* index to start on */
    unsigned int innerIndex = 0;
    for (unsigned int i = 0; i < context.threadCount_; ++i) {

        /* each thread has its values set*/
        threadContext *args = malloc(sizeof(struct threadContext));
        args->start = innerIndex; // what index in the inner array (excluding border array) should we start from, 3 x 3 = 0 to 8 index
        args->blockSize = context.operationsPerThread[i]; // number of operations for that thread to complete
        args->context = &context;
        args->id = i;

        pthread_create(&context.threads[i], 0, (void *(*)(void *)) work, (void *) args);

        // block size is added so the next thread knows where to start from
        innerIndex += context.operationsPerThread[i];
    }

    // main thread controls the sequential part of the threads that have started running
    controlThread(&context);

    for (unsigned int i = 0; i < context.threadCount_; ++i) {
        pthread_join(context.threads[i], 0);
    }
    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / (double) 1000000000.0;

    print_array(outputBuffer, context.arraySize_);
    printf("elapsed time: %f\n", elapsed);
    return 0;
}