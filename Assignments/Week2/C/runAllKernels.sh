#
echo "Running 4x? Kernels"
make JI_4x1Kernel
make JI_4x2Kernel
make JI_4x4Kernel
make JI_4x8Kernel
make JI_4x12Kernel

echo "Running 8x? Kernels"
make JI_8x1Kernel
make JI_8x2Kernel
make JI_8x4Kernel
make JI_8x6Kernel
make JI_8x8Kernel

echo "Running 12x? Kernels"
make JI_12x1Kernel
make JI_12x2Kernel
make JI_12x4Kernel
make JI_12x6Kernel

echo "Done"
pause
