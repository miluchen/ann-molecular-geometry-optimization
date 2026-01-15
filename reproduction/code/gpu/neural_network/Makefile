train:
	nvcc -arch=sm_37 -dc main.cu cuchemnet.cu trainer.cu dataset.cu config.cu -std=c++11
	nvcc -arch=sm_37 main.o cuchemnet.o trainer.o dataset.o config.o -lcudadevrt -lcublas -o neural_train
train_grad:
	nvcc -arch=sm_37 -dc main_grad.cu cuchemnet.cu trainer.cu dataset.cu config.cu -std=c++11
	nvcc -arch=sm_37 main_grad.o cuchemnet.o trainer.o dataset.o config.o -lcudadevrt -lcublas -o neural_train_grad
test:
	nvcc -arch=sm_37 -dc test.cu cuchemnet.cu trainer.cu dataset.cu config.cu -std=c++11
	nvcc -arch=sm_37 test.o cuchemnet.o trainer.o dataset.o config.o -lcudadevrt -lcublas -o neural_test
