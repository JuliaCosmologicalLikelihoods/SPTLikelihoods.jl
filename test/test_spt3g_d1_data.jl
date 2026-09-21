@testset "SPT-3G D1 fixed data" begin
    ells = collect(2:5)
    windows = [reshape(fill(1.0 / (index + 1), length(ells) * (index + 1)), length(ells), index + 1)
               for index in 1:2]
    data = collect(1.0:5.0)
    covariance = Matrix{Float64}(I, length(data), length(data))

    d = SPT3GD1Data(
        data,
        covariance,
        windows,
        ells;
        spectrum_order=["TT 90x90", "TE 90x90"],
        bins_per_spectrum=[2, 3],
    )

    @test d.data_vector == data
    @test d.ells == ells
    @test d.bins_per_spectrum == [2, 3]
    @test size(d.cov_chol) == (5, 5)
    @test all(size(window, 1) == length(ells) for window in d.windows)

    # Validation tests
    @test_throws DimensionMismatch SPT3GD1Data(
        data, covariance, windows, ells;
        spectrum_order=["TT 90x90", "TE 90x90"], bins_per_spectrum=[1, 3],
    )
    @test_throws DimensionMismatch SPT3GD1Data(
        data, covariance, windows[1:1], ells;
        spectrum_order=["TT 90x90", "TE 90x90"], bins_per_spectrum=[2, 3],
    )
    @test_throws DimensionMismatch SPT3GD1Data(
        data[1:4], covariance, windows, ells;
        spectrum_order=["TT 90x90", "TE 90x90"], bins_per_spectrum=[2, 3],
    )
    @test_throws ArgumentError SPT3GD1Data(
        data, covariance, windows, [2, 4, 5, 6];
        spectrum_order=["TT 90x90", "TE 90x90"], bins_per_spectrum=[2, 3],
    )
    # Non-finite data
    bad_data = copy(data); bad_data[1] = NaN
    @test_throws ArgumentError SPT3GD1Data(bad_data, covariance, windows, ells; spectrum_order=["TT 90x90", "TE 90x90"], bins_per_spectrum=[2, 3])
    # Non-finite covariance
    bad_cov = copy(covariance); bad_cov[1, 1] = Inf
    @test_throws ArgumentError SPT3GD1Data(data, bad_cov, windows, ells; spectrum_order=["TT 90x90", "TE 90x90"], bins_per_spectrum=[2, 3])
    # Non-symmetric covariance
    asym_cov = copy(covariance); asym_cov[1, 2] = 0.5
    @test_throws ArgumentError SPT3GD1Data(data, asym_cov, windows, ells; spectrum_order=["TT 90x90", "TE 90x90"], bins_per_spectrum=[2, 3])
    # Non-positive definite covariance
    nonpos_cov = Matrix{Float64}(-I, length(data), length(data))
    @test_throws ArgumentError SPT3GD1Data(data, nonpos_cov, windows, ells; spectrum_order=["TT 90x90", "TE 90x90"], bins_per_spectrum=[2, 3])
end
