load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")
load("@bazel_tools//tools/build_defs/repo:git.bzl", "git_repository")
load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")

http_archive(
    name = "rules_foreign_cc",
    # TODO: Get the latest sha256 value from a bazel debug message or the latest
    #       release on the releases page: https://github.com/bazelbuild/rules_foreign_cc/releases
    #
    # sha256 = "...",
    strip_prefix = "rules_foreign_cc-f01fd353ee2adcd55aab899c12fa2733223228a1",
    url = "https://github.com/bazelbuild/rules_foreign_cc/archive/f01fd353ee2adcd55aab899c12fa2733223228a1.tar.gz",
)

load("@rules_foreign_cc//foreign_cc:repositories.bzl", "rules_foreign_cc_dependencies")

rules_foreign_cc_dependencies()

_ALL_CONTENT = """\
filegroup(
    name = "all_srcs",
    srcs = glob(["**"]),
    visibility = ["//visibility:public"],
)
"""

# soplex source code repository
http_archive(
    name = "soplex",
    build_file_content = _ALL_CONTENT,
    strip_prefix = "soplex-master",
    urls = [
        "https://github.com/scipopt/soplex/archive/refs/heads/master.zip",
    ],
    # sha256 = "0b8e7465dc5e98c757cc3650a20a7843ee4c3edf50aaf60bb33fd879690d2c73",
)

# Extra Bazel utilities that may be handy.
git_repository(
    name = "bazel_skylib",
    commit = "2ec2e6d",  # release 1.0.3
    remote = "https://github.com/bazelbuild/bazel-skylib.git",
)

# Support for building C++ targets in Bazel.
http_archive(
  name = "rules_cc",
  urls = ["https://github.com/bazelbuild/rules_cc/archive/40548a2974f1aea06215272d9c2b47a14a24e556.zip"],
  strip_prefix = "rules_cc-40548a2974f1aea06215272d9c2b47a14a24e556",
)

# Support for building Python targets in Bazel.
git_repository(
    name = "rules_python",
    commit = "c8c79aa",  # 0.1.0
    remote = "https://github.com/bazelbuild/rules_python.git",
)

# Google unit testing framework.
# See: https://google.github.io/googletest/primer.html for introduction.
http_archive(
  name = "com_google_googletest",
  urls = ["https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip"],
  strip_prefix = "googletest-609281088cfefc76f9d0ce82e1ff6c30cc3591e5",
)

# Google protocol buffers (structs on "steroids", useful to serialize to files
# and on wire).
git_repository(
    name = "com_google_protobuf",
    commit = "983d115",  # release v3.15.3
    remote = "https://github.com/protocolbuffers/protobuf.git",
)
load("@com_google_protobuf//:protobuf_deps.bzl", "protobuf_deps")

# Google abseil library with many useful utilities (e.g., flags, status).
# See: https://abseil.io/docs/cpp/guides/
#
# Note, some utilities are not yet (fully) covered (e.g., logging) and are thus
# provided in OR-Tools below.
git_repository(
    name = "com_google_absl",
    commit = "278e0a0",  # release 20210324.2
    remote = "https://github.com/abseil/abseil-cpp.git",
)

# Google OR-Tools. Contains many LP/MIP related utilities, including GLOP
# (LP solver) and SCIP (MIP solver).
#
# Note, OR-Tools includes some of the Google utilities not yet (fully) included
# in abseil (https://github.com/google/or-tools/tree/stable/ortools/base).
# In particular, OR-Tools includes logging routines, which are not part of
# abseil yet. Importantly, we do not depend on and use "glog", because it's old
# and depends on  "gflags" (now replaced by abseil/Flags).
git_repository(
    name = "com_google_ortools",
    commit = "fa84bc0",  # release v8.2
    remote = "https://github.com/google/or-tools.git",
)

# zlib -- needed to compile SCIP.
http_archive(
    name = "zlib",
    build_file = "@com_google_protobuf//:third_party/zlib.BUILD",
    sha256 = "c3e5e9fdd5004dcb542feda5ee4f0ff0744628baf8ed2dd5d66f8ca1197cb1a1",
    strip_prefix = "zlib-1.2.11",
    urls = [
        "https://mirror.bazel.build/zlib.net/zlib-1.2.11.tar.gz",
        "https://zlib.net/zlib-1.2.11.tar.gz",
    ],
)

# bliss -- needed to compile SCIP with symmetry breaking.
http_archive(
    name = "bliss",
    build_file = "@com_google_ortools//bazel:bliss.BUILD",
    patches = ["@com_google_ortools//bazel:bliss-0.73.patch"],
    sha256 = "f57bf32804140cad58b1240b804e0dbd68f7e6bf67eba8e0c0fa3a62fd7f0f84",
    url = "http://www.tcs.hut.fi/Software/bliss/bliss-0.73.zip",
)

# SCIP -- MIP solver, available via OR-Tools MPSolver. Also includes SoPlex. (It does not)
http_archive(
    name = "scip",
    build_file = "@com_google_ortools//bazel:scip.BUILD",
    patches = ["@com_google_ortools//bazel:scip.patch"],
    sha256 = "033bf240298d3a1c92e8ddb7b452190e0af15df2dad7d24d0572f10ae8eec5aa",
    url = "https://github.com/google/or-tools/releases/download/v7.7/scip-7.0.1.tgz",
)
