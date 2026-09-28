using TestItemRunner

# Select test items by tag: `Pkg.test(test_args=["joint_diagonalization"])` runs only items
# tagged :joint_diagonalization, "noaqua" skips items tagged :aqua, and "all" also runs the
# :slow items, which are skipped by default. Arguments can also be passed as
# ORBITALROTATIONSTK_TEST_ARGS="a-b-c".
args = isempty(ARGS) ?
    split(get(ENV, "ORBITALROTATIONSTK_TEST_ARGS", ""), "-"; keepempty=false) : ARGS
included = [Symbol(arg) for arg in args if arg != "all" && !startswith(arg, "no")]
excluded = [Symbol(arg[3:end]) for arg in args if startswith(arg, "no")]
# :slow items run on request only: via "all", or when tags are selected explicitly
("all" in args || !isempty(included)) || push!(excluded, :slow)

function testfilter(ti)
    any(in(ti.tags), excluded) && return false
    return isempty(included) || any(in(ti.tags), included)
end

println("Running OrbitalRotationsTK tests")
isempty(included) || println("    Included tags: ", join(included, ", "))
isempty(excluded) || println("    Excluded tags: ", join(excluded, ", "))

@run_package_tests filter=testfilter verbose=true
