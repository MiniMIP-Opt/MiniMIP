# MiniMIP

MiniMIP is an open source, machine learning oriented Mixed-Integer Programming (MIP) solver.
We provide a range of interfaces for all aspects of solving MIPs (e.g. heuristics, cut generators, LP solvers), 
supplying users with a constant view of the internal state and allowing them to propose modifications that are integrated into the global state internally.

This project is still under heavy development and is not ready for practical application.

# Build system

## Bazel
The `bazel` build environment composes of:
* The file `WORKSPACE` that tells from where to fetch all the dependencies
* The file `.bazelrc` that defines commonly used options while building
(e.g., in MiniMIP we always compile with a flag `--cxxopt=-std=c++2a`)
* The build rules for all targets, which are stated explicitly in each
(sub-)directory's `BUILD.bazel` file.

If you don't have Bazel installed already, you can download it from https://bazel.build/.

To build a single target in `minimip/examples` run:
```
bazel build //minimip/examples:example
```

To build all targets in `minimip/examples` run:

```
bazel build //minimip/examples:all
```

Bazel puts all artifacts generated during compilation in the top directory in
directories that start with `bazel-`. The most useful is `bazel-bin` -- this is
where the executables are kept (including the executables for unit tests).

Since you have built your binaries already, you can now run the test:
```
bazel-bin/minimip/examples/example_test
```

Or the main executable:
```
bazel-bin/minimip/examples/example_main --logtostderr
```

The `--logtostderr` flag tells the logging system to output logs to `stderr`
(by default is streamed to a file).

You could achieve the same thing by running from bazel directly:
```
bazel run //minimip/examples:example_main -- -logtostderr
```
The `-- -logtostderr` is needed to pass the run-time argument to `example_main`.
All command-line arguments after "--" are passed to the launched binary.

Similarly, you can also run the tests directly from bazel:
```
bazel test //minimip/examples:example_test
```

Or, if you have many tests to run (we don't have here, so the effect will be
identical as above):
```
bazel test //minimip/examples:all
```

Bazel by default builds in "fast-build" mode. To maximize binary performance
(e.g., for benchmarking purposes) compile in "opt" mode:
```
bazel build -c opt //minimip/example:example_main
```
To include symbols allowing for debugging (with `gdb`) compile in `dbg` mode:
```
bazel build -c dbg //minimip/example:example_main
```
See the [bazel doc](https://docs.bazel.build/versions/master/user-manual.html#flag--compilation_mode)
for more details.

## CMake
We also provide a `CMakeLists.txt` file for building the project with CMake.
To build the project with CMake, you can run the following commands:
```
cmake -S. -Bbuild
cmake --build build -v -j 12
(cd build && ctest)
```

## Dependencies
See `WORKSPACE` or `CMakeLists.txt` file for a list of dependencies and extra information. These
dependencies are automatically pulled by `bazel`.

**Important**
Boost still needs to be installed on the machine for SoPlex to read LP files and be tested.
On Debian machines, the package `libboost-all-dev` covers it.

# Style guide
Google C++ style guide is available at:
https://google.github.io/styleguide/cppguide.html

To automatically format the syntax you can run:
`clang-format -i PATH_TO_FILES/*.h PATH_TO_FILES/*.cc` which will pick up the style specified in `.clang-format`.

Be careful not run clang-format on BUILD files. To format bazel's BUILD files
automatically, run `buildifier` (see https://github.com/bazelbuild/buildtools).

# VSCode setup
There is an official VSCode plugin for bazel called `vscode-bazel`. While it allows context highlights when editing Bazel files (`BUILD`, `WORKSPACE`, etc.) it doesn't configure IntelliSense for code files compiled using Bazel. Instead, use the following steps to get code suggestions using `clangd`. Tested on Linux Mint 20 (Ubuntu 20.04 derivative) on 2022-10-04.

1. Make sure `clangd` is installed on your system
```
sudo apt install clangd
```
2. Download the [Bazel Compilation Database](https://github.com/grailbio/bazel-compilation-database) tool.
3. Go to MiniMIP's root directory (where the `WORKSPACE` file is located) and run the `generate.py` script you just downloaded. When the script has completed (may take several minutes), you should see a `compile_commands.json` file.
4. Install the `vscode-clangd` extension for VSCode.

    Note: The clangd extension is incompatible with Microsoft's IntelliSense, and will ask you to disable it. If you want to control this manually, the setting is found under C/C++: Intelli Sense Engine.
5. In your VSCode settings, find Clangd: Arguments and add `--compile-commands-dir=${workspaceFolder}/` and make sure Clangd: Check Updates is enabled.
6. Restart VSCode.

Your MiniMIP files should now give you contextual highlighting, symbol lookup, clang-tidy warnings and more. You may need to occasionally repeat step 3 above if the information in `compile_commands.json` becomes stale, e.g. if you add a new target or introduce new dependencies.
