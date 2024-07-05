load("@rules_cc//cc:defs.bzl", "cc_library", "cc_binary")

cc_library(
    name = "libsoplex",
    srcs = glob(
      [
        "src/soplex/*.cpp",
        "src/soplex/external/fmt/*.cc",
      ],
    ),
    copts = [
        "-w",
        "-fexceptions",
        "-DSOPLEX_NO_CONFIG_HEADER",
    ],
    defines = ["SOPLEX_NO_CONFIG_HEADER"],
    includes = [
    "src",
    "src/soplex/external",
    ],
    textual_hdrs = glob([
        "src/soplex/*.h",
        "src/soplex/*.hpp",
        "src/soplex/external/fmt/*.h",
    ]) + [
        "src/soplex.h",
        "src/soplex.hpp",
        "src/soplex/git_hash.cpp",
    ],
    visibility = ["//visibility:public"],
    deps = [],
)

cc_binary(
    name = "soplex",
    srcs = [
        "src/soplexmain.cpp",
    ],
    copts = [
        "-w",
        "-fexceptions",
    ],
    deps = [
        ":libsoplex",
    ],
)
