@testset "SPT-3G D1 nuisance parameter container" begin
    fiducial = SPT3GD1Parameters()

    # 1. The vector boundary has the released width and is fully finite
    v = Vector(fiducial)
    @test length(v) == 43
    @test length(SPT3G_D1_PARAMETER_NAMES) == 43
    @test all(isfinite, v)
    @test v == collect(fiducial)
    @test eltype(v) === Float64

    # 2. Vector(p) is ordered exactly as SPT3G_D1_PARAMETER_NAMES, which pins every
    #    field name and slot against the canonical release ordering.
    @test v == [Float64(SPT3G_D1_FIDUCIAL_PARAMETERS[Symbol(name)])
                for name in SPT3G_D1_PARAMETER_NAMES]

    # 3. Vector <-> SPT3GD1Parameters round-trips. A strictly increasing probe makes
    #    any duplicated, swapped, or dropped slot observable.
    probe = collect(1.0:43.0)
    @test Vector(SPT3GD1Parameters(probe)) == probe

    # 4. The NamedTuple constructor agrees with the zero-argument one
    @test Vector(SPT3GD1Parameters(SPT3G_D1_FIDUCIAL_PARAMETERS)) == v

    # 5. Every parameter actually reaches a named field rather than a silent default
    perturbed = SPT3GD1Parameters(v .+ 0.5)
    @test count(!=(0.0), Vector(perturbed) .- v) == 43

    # 6. Width is enforced at the vector boundary
    @test_throws DimensionMismatch SPT3GD1Parameters(collect(1.0:42.0))
    @test_throws DimensionMismatch SPT3GD1Parameters(collect(1.0:44.0))

    # 7. The container is differentiable-parameter generic, not Float64-locked
    dual_vector = ForwardDiff.Dual.(v, 1.0)
    @test Vector(SPT3GD1Parameters(dual_vector)) == dual_vector
end
