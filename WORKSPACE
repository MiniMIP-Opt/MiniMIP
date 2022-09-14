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

# Extra Bazel utilities that may be handy.
http_archive(
    name = "bazel_skylib",
    strip_prefix = "bazel-skylib-1.0.3",
    urls = [
        "https://github.com/bazelbuild/bazel-skylib/archive/refs/tags/1.0.3.zip",
    ],
    sha256 = "28f81e36692e1d87823623a99966b2daf85af3fdc1b40f98e37bd5294f3dd185",
)

# Support for building C++ targets in Bazel.
http_archive(
  name = "rules_cc",
  urls = ["https://github.com/bazelbuild/rules_cc/archive/40548a2974f1aea06215272d9c2b47a14a24e556.zip"],
  strip_prefix = "rules_cc-40548a2974f1aea06215272d9c2b47a14a24e556",
)

# Support for building Python targets in Bazel.
http_archive(
    name = "rules_python",
    urls = [
        "https://github.com/bazelbuild/rules_python/archive/refs/tags/0.1.0.zip",
    ],
    strip_prefix = "rules_python-0.1.0",
    sha256 = "e9afc70f3d14151f1f8759dc224159b029068a66d11f554e0df814d001a345d9",
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
http_archive(
    name = "com_google_protobuf",
    urls = [
        "https://github.com/protocolbuffers/protobuf/archive/refs/tags/v3.19.4.zip",
    ],
    strip_prefix = "protobuf-3.19.4",
    sha256 = "25680843adf0c3302648d35f744e38cc3b6b05a6c77a927de5aea3e1c2e36106",
)
# Load common dependencies.
load("@com_google_protobuf//:protobuf_deps.bzl", "protobuf_deps")
protobuf_deps()

# Google abseil library with many useful utilities (e.g., flags, status).
# See: https://abseil.io/docs/cpp/guides/
#
# Note, some utilities are not yet (fully) covered (e.g., logging) and are thus
# provided in OR-Tools below.
http_archive(
    name = "com_google_absl",
    urls = ["https://github.com/abseil/abseil-cpp/archive/refs/tags/20220623.1.zip"],
    strip_prefix = "abseil-cpp-20220623.1",
    sha256 = "54707f411cb62a26a776dad5fd60829098c181700edcd022ea5c2ca49e9b7ef1",
)

# Abseil uses this to check the current platform of the system. It must be
# explicitly loaded here, since Jenkins doesn't do it automatically for some
# reason.
http_archive(
    name = "platforms",
    sha256 = "a879ea428c6d56ab0ec18224f976515948822451473a80d06c2e50af0bbe5121",
    strip_prefix = "platforms-da5541f26b7de1dc8e04c075c99df5351742a4a2",
    urls = ["https://github.com/bazelbuild/platforms/archive/da5541f26b7de1dc8e04c075c99df5351742a4a2.zip"],  # 2022-05-27
)

# Google OR-Tools. Contains many LP/MIP related utilities, including GLOP
# (LP solver) and SCIP (MIP solver).
#
# Note, OR-Tools includes some of the Google utilities not yet (fully) included
# in abseil (https://github.com/google/or-tools/tree/stable/ortools/base).
# In particular, OR-Tools includes logging routines, which are not part of
# abseil yet. Importantly, we do not depend on and use "glog", because it's old
# and depends on  "gflags" (now replaced by abseil/Flags).
http_archive(
    name = "com_google_ortools",
    urls = ["https://github.com/google/or-tools/archive/refs/tags/v9.4.zip"],
    strip_prefix = "or-tools-9.4",
    sha256 = "c547ac48ac2605e42ec9b003a0adf9345fbe04831fa0858b0d7479125480f93a",
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

http_archive(
    name = "scip",
    build_file = "@com_google_ortools//bazel:scip.BUILD",
    patches = ["@com_google_ortools//bazel:scip.patch"],
    patch_args = ["-p1"],
    sha256 = "ed5535c5def3ebb29cf12ae0309c88e679a6bfee97f7b314b3c342fc7dfbf083",
    strip_prefix = "scip-801",
    url = "https://github.com/scipopt/scip/archive/refs/tags/v801.zip",
)

# Eigen is needed by ortools to compile the PDLP solver, which isn't used by
# MiniMIP but is included in targets we depend on.
http_archive(
    name = "eigen",
    url = "https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.zip",
    strip_prefix = "eigen-3.4.0",
    sha256 = "1ccaabbfe870f60af3d6a519c53e09f3dcf630207321dffa553564a8e75c4fc8",
    build_file_content =
"""
cc_library(
    name = 'eigen3',
    srcs = [],
    includes = ['.'],
    hdrs = glob(['Eigen/**']),
    visibility = ['//visibility:public'],
)
""")

# SoPlex -- available from source code repository and compiled via `cmake` rule.
# TODO(lpawel): See if we can get SoPlex from OR-Tools (since it contains SCIP).
_ALL_CONTENT = """\
filegroup(
    name = "all_srcs",
    srcs = glob(["**"]),
    visibility = ["//visibility:public"],
)
"""

http_archive(
    name = "soplex",
    build_file_content = _ALL_CONTENT,
    strip_prefix = "soplex-master",
    urls = [
        "https://github.com/scipopt/soplex/archive/refs/heads/master.zip",
    ],
    # sha256 = "0b8e7465dc5e98c757cc3650a20a7843ee4c3edf50aaf60bb33fd879690d2c73",
)
