#pragma offload_attribute(push, target(mic))
#include <sstream>
#include <string>
#include <iostream>
using namespace std;

#include <unistd.h>
#include <sys/syscall.h>
#include <pthread.h>
pid_t gettid(void)
{
    return syscall(SYS_gettid);
}
#pragma offload_attribute(pop)

extern "C"
{
    __declspec(target(mic))
    void print_affinity(int team_id, int thread_id)
    {
        cpu_set_t cpuset;
        pthread_t thread;
        thread = pthread_self();
        CPU_ZERO(&cpuset);
        pthread_getaffinity_np(thread, sizeof(cpu_set_t), &cpuset);
        ostringstream os;
        os << "Team " << team_id << ", Thread " << thread_id << ": { ";
        for (int c = 0; c < CPU_SETSIZE; c++)
        {
            if (CPU_ISSET(c, &cpuset))
            {
                os << c << " ";
            }
        }
        os << "}";
        cout << os.str() << endl;
    }

}
