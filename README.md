### MPI Project (`mpi.c`)

The MPI project utilizes the MPI (Message Passing Interface) library for parallel computing, enabling processes to communicate with each other in a distributed computing environment. The code structure demonstrates the implementation of array relaxation using MPI for parallel processing. Key aspects include:

- Initialization of the MPI environment and handling errors during MPI initialization.
- Definition of a `context_t` structure to encapsulate the environment setup, including array size, precision, rank of processors, and block sizes.
- Implementation of array relaxation algorithm, which includes scattering the input buffer across processors, performing computation to relax array elements, and gathering the results back.
- Use of MPI functions such as `MPI_Init`, `MPI_Comm_rank`, `MPI_Comm_size`, `MPI_Scatterv`, `MPI_Gatherv`, `MPI_Bcast`, `MPI_Reduce`, and `MPI_Finalize` to manage parallel tasks and communication.
- Calculation of elapsed time using `MPI_Wtime` to measure the performance of the parallel processing.

### Pthreads Project (`pthread.c`)

The Pthreads project employs POSIX threads (Pthreads) to achieve parallel computing by dividing the work among multiple threads running concurrently within a single process. The implementation showcases a parallel computation technique using threads to perform operations on a matrix. Highlights include:

- Creation of a `Context` structure to store thread-related data such as thread count, array size, precision, and buffers for input and output.
- Implementation of a function to print array values for debugging and visualization purposes.
- Calculation of averages using neighboring cell values to demonstrate data dependency handling between threads.
- Management of thread synchronization using `pthread_barrier_wait` to coordinate the execution flow among multiple threads.
- Division of work among threads, each responsible for a segment of the array, and aggregation of results after computation.
- Measurement of elapsed time using `clock_gettime` to evaluate the performance benefits of multithreading.
