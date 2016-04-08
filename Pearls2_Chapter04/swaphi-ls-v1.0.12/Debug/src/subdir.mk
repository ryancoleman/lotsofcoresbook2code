################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Align.cpp \
../src/AlignCore.cpp \
../src/AlignCoreNaive.cpp \
../src/AlignCoreTiling.cpp \
../src/SeqFileParser.cpp \
../src/Sequence.cpp \
../src/Utils.cpp \
../src/main.cpp 

OBJS += \
./src/Align.o \
./src/AlignCore.o \
./src/AlignCoreNaive.o \
./src/AlignCoreTiling.o \
./src/SeqFileParser.o \
./src/Sequence.o \
./src/Utils.o \
./src/main.o 

CPP_DEPS += \
./src/Align.d \
./src/AlignCore.d \
./src/AlignCoreNaive.d \
./src/AlignCoreTiling.d \
./src/SeqFileParser.d \
./src/Sequence.d \
./src/Utils.d \
./src/main.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	mpic++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


