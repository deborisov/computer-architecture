//Borisov Dmitry csbse192 
//Variant 4
#include <iostream>
#include <thread>
#include <string>
#include <fstream>
#include <iomanip>
#include <pthread.h>
#include <semaphore.h>
#include <random>

enum Status {
	THINKING = 0,
	HUNGRY = 1,
	EATING = 2
};

const int N = 5;

Status state[N];
int philosophers[N];
sem_t mutex;
sem_t display;
sem_t phil_sem[N];

void think(int phil_num, int time) {
	sem_wait(&display);
	std::cout << "Philosopher " << phil_num + 1 << " is thinking..." << std::endl;
	sem_post(&display);
	std::this_thread::sleep_for(std::chrono::milliseconds(time));
}

void eat(int phil_num, int time) {
	sem_wait(&display);
	std::cout << "Philosopher " << phil_num + 1 << " started his meal :)" << std::endl;
	sem_post(&display);
	std::this_thread::sleep_for(std::chrono::milliseconds(time));
	sem_wait(&display);
	std::cout << "Philosopher " << phil_num + 1 << " ended his meal :)" << std::endl;
	sem_post(&display);
}

void try_begin_to_eat(int phil_num) {
	if (state[phil_num] == HUNGRY &&
		state[(phil_num + N - 1) % N] != EATING &&
		state[(phil_num + 1) % N] != EATING) {
		state[phil_num] = EATING;
		sem_post(&phil_sem[phil_num]);
	}
}

void take_forks(int phil_num) {
	sem_wait(&mutex);
	state[phil_num] = HUNGRY;
	sem_wait(&display);
	std::cout << "Philosopher " << phil_num + 1 << " is hungry :(" << std::endl;
	sem_post(&display);
	try_begin_to_eat(phil_num);
	sem_post(&mutex);
	sem_wait(&phil_sem[phil_num]);
}

void put_forks(int phil_num) {
	sem_wait(&mutex);
	state[phil_num] = THINKING;
	try_begin_to_eat((phil_num + N - 1) % N);
	try_begin_to_eat((phil_num + 1) % N);
	sem_post(&mutex);
}

void* philosopher(void* number) {
	int phil_num = *((int*)number);
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> dis(1000, 5000);
	while (true) {
		think(phil_num, dis(gen));
		take_forks(phil_num);
		eat(phil_num, dis(gen));
		put_forks(phil_num);

	}
}

int main(int argc, char* argv[])
{
	if (argc != 2) {
		std::cout << "Please input cmd argument - time to execute in seconds.";
		return 1;
	}
	int time_to_sleep = std::stoi(argv[1]);
	for (int i = 0; i < N; ++i) {
		philosophers[i] = i;
	}
	int i;
	pthread_t thread_id[N];
	sem_init(&mutex, 0, 1);
	sem_init(&display, 0, 1);
	for (i = 0; i < N; i++)
		sem_init(&phil_sem[i], 0, 0);
	for (i = 0; i < N; i++) {
		pthread_create(&thread_id[i], nullptr,
			philosopher, &philosophers[i]);
	}
	std::this_thread::sleep_for(std::chrono::seconds(time_to_sleep));
}

