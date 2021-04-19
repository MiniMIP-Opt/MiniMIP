# MiniMIP project

A modular, minimalistic, robust MILP solver.

# Build system

## Bazel
The `bazel` build environment composes of:
* The file `WORKSPACE` that tells from where to fetch all the dependencies
* The gile `.bazelrc` that defines commonly used options while building (e.g.,
in MiniMIP we alsays complile with a flag `--cxxopt=-std=c++17`
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

Similary, you can also run the tests directly from bazel:
```
bazel test //minimip/examples:example_test
```

Or, if you have many tests to run (we don't have here, so the effect will be
identical as above):
```
bazel test //minimip/examples:all
```

Bazel by default builds in "fastbuild" mode. To maximize binary performance
(e.g., for benchmarking purposes) compile in "opt" mode:
```
bazel build -c opt //minimip/example:example_main
```
To include symbols allowing for debugging (with `gdb`) complile in `dbg` mode:
```
bazel build -c dbg //minimip/example:example_main
```
See [bazel doc](
https://docs.bazel.build/versions/master/user-manual.html#flag--compilation_mode)
for more details.


## Dependencies
See WORKSPACE file for a list of dependencies and extra information. These
dependencies are automatically pulled by `bazel`.

# Style guide
Google C++ style guide is available at:
https://google.github.io/styleguide/cppguide.html

To automatically format the syntax you can run
```clang-format --style=google -i PATH_TO_FILES```
on all files you edit / add. This will reformat the file in-place (due to `-i`)
to conform to Google coding style (due to `--style=google`).

Be careful not run clang-format on BUILD files. To format bazel's BUILD files
automatically, run `buildifier` (see https://github.com/bazelbuild/buildtools).

