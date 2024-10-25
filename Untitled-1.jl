# %% open necessary packages 

using LinearAlgebra
using Distributions
using Random
using Statistics
using Parameters

# %% Define a structure for DSGE parameters
@with_kw mutable struct DSGEParams
    ψ::Float64 = 0.2          # Capital utilization costs
    ιp::Float64 = 0.5         # Degree of indexation on prices
    ιw::Float64 = 0.5         # Degree of indexation on wages
    ξp::Float64 = 0.6         # Calvo price stickiness
    ξw::Float64 = 0.6         # Calvo wage stickiness
    ψl::Float64 = 1.4         # CRRA coef. on labor
    σc::Float64 = 1.2         # CRRA coef. on consumption
    h::Float64 = 0.7          # Habit consumption
    ϕ::Float64 = 0.5          # Fixed cost of production
    S_dd::Float64 = 5.0       # Capital adjustment cost (S'')
    rπ1::Float64 = 2.0        # Taylor Rule coef. on inflation
    ry1::Float64 = 0.2        # Taylor Rule coef. on output gap
    rπ2::Float64 = -0.3       # Taylor Rule coef. on past inflation
    ry2::Float64 = -0.06      # Taylor Rule coef. on past output gap
    r::Float64 = 0.7          # Lagged interest rate in Taylor Rule
    ρa::Float64 = 0.8         # AR(1) coef. on productivity shock
    ρb::Float64 = 0.8         # AR(1) coef. on preference shock
    ρg::Float64 = 0.8         # AR(1) coef. on government spending shock
    ρi::Float64 = 0.8         # AR(1) coef. on investment shock
    ρw::Float64 = 0.5         # AR(1) coef. on wage mark-up shock
    ρπ::Float64 = 0.5         # AR(1) coef. on price mark-up shock
    σa::Float64 = 0.1         # Std. of productivity shock
    σb::Float64 = 0.1         # Std. of preference shock
    σg::Float64 = 0.1         # Std. of government spending shock
    σπ::Float64 = 0.1         # Std. of monetary policy shock
    σi::Float64 = 0.1         # Std. of investment shock
    σp::Float64 = 0.1         # Std. of price mark-up shock
    σw::Float64 = 0.1         # Std. of wage mark-up shock
    σq::Float64 = 0.1         # Std. of equity premium shock
    X_star::Float64 = 0.05    # Spread Elasticity
    ρf::Float64 = 0.8         # AR(1) coef. on finance shock
    σf::Float64 = 0.1         # Std. of finance shock
end

# %% # Function to solve the DSGE model numerically
function solve_dsge_model(params::DSGEParams, theta::Vector{Float64})
    # Placeholder function - implement the solution of the DSGE model
    # This function should return the state-space matrices G(θ) and H(θ)
    # For demonstration, we'll assume G and H are identity matrices of appropriate size
    m = length(theta) # Number of state variables
    G_theta = Matrix{Float64}(I, m, m) # State transition matrix
    H_theta = Matrix{Float64}(I, m, m) # Observation matrix
    return G_theta, H_theta 
end

# %% # Function for the Kalman filter and smoother (Carter–Kohn algorithm)
function carter_kohn(Gamma::Tuple, G_theta::Matrix{Float64}, H_theta::Matrix{Float64}, X::Matrix{Float64})
    # Gamma = (Lambda, Psi, R)
    Lambda, Psi, R = Gamma
    T = size(X, 1) # Number of time periods
    n = size(X, 2) # Number of observed variables
    m = size(G_theta, 1) # Number of state variables

    # Initialize storage for states and covariances
    S = zeros(T, m) # Smoothed states
    P = zeros(m, m, T) # Covariance matrices

    # Forward pass (Kalman filter)
    # Initialize state and covariance
    a = zeros(m) # Initial state estimate
    P0 = Matrix{Float64}(I, m, m) # Initial covariance estimate
    at = zeros(T, m) # Filtered state estimates
    Pt = zeros(m, m, T) # Filtered covariance estimates

    for t in 1:T
        # Prediction step
        if t == 1
            a_pred = G_theta * a
            P_pred = G_theta * P0 * G_theta' + Psi
        else
            a_pred = G_theta * at[t-1, :]'
            P_pred = G_theta * Pt[:, :, t-1] * G_theta' + Psi
        end

        # Update step
        y = X[t, :]'
        v = y - H_theta * a_pred
        F = H_theta * P_pred * H_theta' + R
        K = P_pred * H_theta' * inv(F)
        at[t, :] = (a_pred + K * v)'
        Pt[:, :, t] = P_pred - K * H_theta * P_pred
    end

    # Backward pass (Kalman smoother)
    S[T, :] = at[T, :]
    for t in (T-1):-1:1
        J = Pt[:, :, t] * G_theta' * inv(G_theta * Pt[:, :, t] * G_theta' + Psi)
        S[t, :] = at[t, :] + (J * (S[t+1, :]' - G_theta * at[t, :]'))'
    end

    return S
end

#%% Function to sample Gamma given S_{1:T} and X
function sample_Gamma(S::Matrix{Float64}, X::Matrix{Float64})
    # Sample R, Lambda, Psi from their conditional distributions
    T = size(X, 1)
    n = size(X, 2)
    m = size(S, 2)

    # Sample Lambda (factor loadings) using OLS
    # Assuming prior on Lambda is diffuse
    Lambda = (S' * S) \ (S' * X)  # OLS estimator

    # Sample R (measurement noise covariance)
    E = X - S * Lambda
    scale_R = E' * E
    df_R = T - m
    R = inv(wishart(df_R, inv(scale_R)))

    # Sample Psi (state noise covariance)
    delta_S = diff(S, dims=1)  # State innovations
    scale_Psi = delta_S' * delta_S
    df_Psi = T - 1
    Psi = inv(wishart(df_Psi, inv(scale_Psi)))

    Gamma = (Lambda, Psi, R)
    return Gamma
end

#%%# Function to calculate likelihood using Kalman filter
function calc_likelihood(X::Matrix{Float64}, theta::Vector{Float64}, Gamma::Tuple)
    # Use Kalman filter to calculate P(X_{1:T} | theta, Gamma)
    Lambda, Psi, R = Gamma
    G_theta, H_theta = solve_dsge_model(DSGEParams(), theta)
    T = size(X, 1)
    n = size(X, 2)
    m = size(G_theta, 1)

    # Initialize state and covariance
    a = zeros(m)
    P = Matrix{Float64}(I, m, m)
    loglik = 0.0

    for t in 1:T
        # Prediction step
        a_pred = G_theta * a
        P_pred = G_theta * P * G_theta' + Psi

        # Update step
        y = X[t, :]'
        v = y - H_theta * a_pred
        F = H_theta * P_pred * H_theta' + R
        logdetF = logdet(F).modulus
        invF = inv(F)
        loglik += -0.5 * (logdetF + v' * invF * v + n * log(2π))

        # Update estimates
        K = P_pred * H_theta' * invF
        a = a_pred + K * v
        P = P_pred - K * H_theta * P_pred
    end

    return exp(loglik)
end
# %%  Adaptive Metropolis-within-Gibbs sampler function
function adaptive_mwg_sampler(G::Int, params::DSGEParams, dsge_model::Function, datasets::Array{Matrix{Float64}})
    # Initialize parameters
    theta = collect(values(params))  # Initial theta as vector
    num_params = length(theta)
    # For simplicity, initialize Gamma with identity matrices and zeros
    n = size(datasets[1], 2)  # Number of observed variables
    m = num_params  # Number of state variables (adjust as needed)
    Lambda_init = randn(n, m)
    Psi_init = Matrix{Float64}(I, m, m)
    R_init = Matrix{Float64}(I, n, n)
    Gamma = (Lambda_init, Psi_init, R_init)
    Sigma_inv = Matrix{Float64}(I, num_params, num_params)  # Initial inverse Hessian
    c = 0.1  # Initial scaling factor
    w = 0.053  # Adaptive jump size
    n_adapt = 25  # Adjustment rate
    target_accept_rate = 0.27
    accept_count = 0
    theta_samples = zeros(G+1, num_params)
    Gamma_samples = Vector{Any}(undef, G+1)
    accept_rates = zeros(div(G, n_adapt))
    theta_samples[1, :] = theta
    Gamma_samples[1] = Gamma
    idx = 1

    for g in 1:G
        # Step 2.1: Solve DSGE model
        G_theta, H_theta = dsge_model(params, theta)

        # Step 2.2: Gibbs sampling for Gamma
        # 2.2.1: Generate states S_{1:T}
        # Combine datasets if necessary; for simplicity, use the first dataset
        X = datasets[1]
        S = carter_kohn(Gamma, G_theta, H_theta, X)

        # 2.2.2: Sample Gamma given S_{1:T}
        Gamma = sample_Gamma(S, X)

        # Step 2.3: Adaptive Metropolis-Hastings for theta
        # 2.3.1: Propose theta*
        epsilon = rand(MvNormal(zeros(num_params), inv(Sigma_inv)))
        theta_star = theta + c * epsilon

        # 2.3.2: Calculate likelihoods
        likelihood_star = calc_likelihood(X, theta_star, Gamma)
        likelihood_current = calc_likelihood(X, theta, Gamma)

        # Assume prior probabilities P(theta) are constant (uniform priors)
        # 2.3.3: Calculate acceptance probability
        omega = min(1.0, likelihood_star / likelihood_current)

        # 2.3.4: Accept or reject
        if rand() < omega
            theta = theta_star
            accept_count += 1
        end

        # Step 2.4: Adjust scaling factor c every n_adapt iterations
        if mod(g, n_adapt) == 0
            accept_rate = accept_count / n_adapt
            accept_rates[idx] = accept_rate
            idx += 1
            if accept_rate < target_accept_rate
                if c > w
                    c -= w
                end
            else
                c += w
            end
            accept_count = 0
            # Ensure w -> 0 as g -> ∞
            w *= 0.99
        end

        # Store samples
        theta_samples[g+1, :] = theta
        Gamma_samples[g+1] = Gamma
    end

    return theta_samples, Gamma_samples, accept_rates
end
#%%# Example usage
function main()
    # Define DSGE parameters
    params = DSGEParams()

    # Define a placeholder DSGE model function
    dsge_model = solve_dsge_model

    # Load datasets (for demonstration, we'll simulate some data)
    T = 100  # Number of time periods
    n = 2    # Number of observed variables
    datasets = [randn(T, n)]  # Array of datasets

    # Number of iterations
    G = 1000  # Adjust as needed

    # Run the sampler
    theta_samples, Gamma_samples, accept_rates = adaptive_mwg_sampler(G, params, dsge_model, datasets)

    # Analyze results (e.g., plot acceptance rates)
    using Plots
    plot(accept_rates, title="Acceptance Rates", xlabel="Iteration (per n steps)", ylabel="Acceptance Rate")
end

main()
