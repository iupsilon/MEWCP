################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../MEWCP_combinatorial_bound.c \
../MEWCP_dsdp.c \
../MEWCP_explicit_enumeration.c \
../MEWCP_tabu.c \
../converter_dsdp.c \
../main.c 

OBJS += \
./MEWCP_combinatorial_bound.o \
./MEWCP_dsdp.o \
./MEWCP_explicit_enumeration.o \
./MEWCP_tabu.o \
./converter_dsdp.o \
./main.o 

C_DEPS += \
./MEWCP_combinatorial_bound.d \
./MEWCP_dsdp.d \
./MEWCP_explicit_enumeration.d \
./MEWCP_tabu.d \
./converter_dsdp.d \
./main.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


