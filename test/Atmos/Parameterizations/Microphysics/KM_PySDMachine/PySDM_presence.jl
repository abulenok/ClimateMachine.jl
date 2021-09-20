module PySDMTest

using Dates: print, isequal

using PyCall

using Test

println("[TEST] PySDM presence test")

pysdm = pyimport("PySDM")

@test pysdm isa PyObject
@test pysdm.__name__ == "PySDM"

module_content = py"dir($pysdm)"

@test "Core" in module_content &&
        "physics" in module_content &&
        "builder" in module_content &&
        "initialisation" in module_content &&
        "backends" in module_content

end        