# CSC417-homeworks

This repository is the collection of notes and homeworks based on the course [CSC417-physics-based animation](https://github.com/dilevin/CSC417-physics-based-animation)

## build

```shell
mkdir build && cd build
cmake ..
make
```

- [mass-spring-1D](./a1-mass-spring-1d/README.md)

```shell
./mass-spring-1d
```

```shell
./mass-spring-1d 'rk'
```

```shell
./mass-spring-1d 'be'
```

```shell
./mass-spring-1d 'se'
```

- [mass-spring-3d](./a2-mass-spring-3d/README.md)

```shell
cmake .. -DCMAKE_BUILD_TYPE=Release
make
```

```shell
./mass-spring-3d
```

- [FEM-3D](./a3-finite-elements-3d/README.md)

```shell
cmake .. -DCMAKE_BUILD_TYPE=Release
make
```

```shell
./finite-elements-3d
```
