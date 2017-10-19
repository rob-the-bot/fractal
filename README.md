# Fractal generation

A simple program in C to generate Newton fractal using arbitrary-precision math library.

# Usage

```bash
make
./fract spec
```

The spec file indicates

- resolution
- range?
- range?
- Zoom level
- Number of frames

> [!NOTE]  
> Testing fractal generation with other polynomials requires recompiling the source file, the committed code is using $z^3+1=0$ but can be easily changed to other polynomials by modifying line 19 at [`fract.c`](fract.c). Since $z^3+1=0$ doesn't have 0 as one of the root, it's very suitable for testing the zoom-in functionality.

## Use it in computing clusters

Here we include an example&mdash;`go.sh`of using it on SLURM.

slurm parameters can be modified (not hard coded in fract.c) at line 22 in `go.sh`: `time mpirun -n ${nodes} ./fract spec`.
`spec` is the name of the specification file passed into fract function, it can be in different names as long as it's in the right format.

Comment: spec names can be given in the format like 20px_15px_zoom_1e_-10_20 which mean width is 20, height is 15, zoom in value 1e-10, 20 frames.
For each testing files:

1. 20px_15px_zoom_1e-10_20 -> testing scalibility (zoom to to 1e-190)
2. 20px_15px_zoom_1e-50_20 -> testing scalibility (zoom to to 1e-950),also compare speed with No.1
3. 1024px_768px_zoom_0.95_820 -> test the speed while using double precision (small problems, speed should be fast)
4. 400px_300px_zoom_0.95_820 -> compare speed with No.3
5. 20px_15px_zoom_0.95_820 -> compare speed with No.4
6. 20px_15px_zoom_0.5_820 -> compare with No.5 which will normally be done in 2 sec. But this one, since mpc is used after certain frame, super slow now

# Visualization

`bash genvis.sh $NUMBER_OF_FRAMES` (e.g. `bash genvis.sh 40`), the `$NUMBER_OF_FRAMES` is required per piazza "Input (genvis) command line argument N the number of frames:"


# Development environment

This code requires

- make
- libmpc-dev
- mpich
