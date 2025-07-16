# Spectral.jl
module Spectral

using Printf

# --- 物理常数 ---
const HBAR = 6.582119569e-16  # eV.s
const HBAR_SQ_OVER_2M0 = 3.80998   # eV.Angstrom^2
const KB_EV_K = 8.617333262e-5   # eV/K

# --- 参数结构体 (替代 Python 的类) ---
struct Spectral1D
    prefix::String
    a_lattice::Float64
    efermi::Float64
    reci_lattice::Float64

    # 构造函数
    function Spectral1D(paras::Dict)
        a_lattice = paras["a_lattice"]
        reci_lattice = 2 * π / a_lattice
        new(paras["prefix"], a_lattice, paras["efermi"], reci_lattice)
    end
end

# --- 辅助函数 ---

function bose_einstein_dist(energy_ev::Float64, temperature_k::Float64)
    if temperature_k <= 1e-9 return 0.0 end
    if energy_ev <= 1e-9 return 1.0 / (exp(1e9) - 1) end # 避免能量为0
    beta_energy = energy_ev / (KB_EV_K * temperature_k)
    if beta_energy > 700 return 0.0 end
    # 处理 exp(x) ≈ 1+x 的情况，避免精度损失
    if abs(exp(beta_energy) - 1.0) < 1e-12
        return 1.0 / beta_energy
    end
    return 1.0 / (exp(beta_energy) - 1.0)
end

function fermi_dirac_dist(energy_ev::Float64, efermi_ev::Float64, temperature_k::Float64)
    if temperature_k <= 1e-9
        if abs(energy_ev - efermi_ev) < 1e-9 return 0.5 end
        return energy_ev < efermi_ev ? 1.0 : 0.0
    end
    beta_energy_mu = (energy_ev - efermi_ev) / (KB_EV_K * temperature_k)
    # exp 函数会自动处理溢出，返回 Inf 或 0.0，无需手动检查
    return 1.0 / (exp(beta_energy_mu) + 1.0)
end

function sigma_FM_integrand_1D(
    spec_params::Spectral1D,
    e_k_bare::Function,          # 裸能带结构函数
    q_scalar::Float64,           # 声子波矢
    k_scalar::Float64, hw_eval_ev::Float64, # 电子参数
    hw_phonon_ev::Float64, temp_k::Float64, efermi_ev::Float64, # 声子和系统参数
    h_delta_ev::Float64          # 展宽
)
    e_k_plus_q = e_k_bare(k_scalar - q_scalar)
    n_phonon = bose_einstein_dist(hw_phonon_ev, temp_k)
    f_electron = fermi_dirac_dist(e_k_plus_q, efermi_ev, temp_k)

    num1 = n_phonon + 1.0 - f_electron
    num2 = n_phonon + f_electron

    # 使用复数直接构造分母
    den1 = complex(hw_eval_ev - e_k_plus_q - hw_phonon_ev, h_delta_ev)
    den2 = complex(hw_eval_ev - e_k_plus_q + hw_phonon_ev, h_delta_ev)

    # Julia 会自动处理除以一个很小的复数的情况
    term1 = num1 / den1
    term2 = num2 / den2

    return term1 + term2
end

function calc_sigma_k_e_FM_1D(
    spec_params::Spectral1D,
    e_k_bare::Function,
    k_scalar_Ainv::Float64, hw_eval_ev::Float64, g0_1D_sq_ev2::Float64, hw_phonon_ev::Float64,
    temp_k::Float64, efermi_ev::Float64, num_q_points::Int,
    h_delta_ev::Float64
)
    q_bz_ws_boundary = 0.5 * spec_params.reci_lattice
    
    # 定义q点网格
    delta_q_for_sum = (2 * q_bz_ws_boundary) / num_q_points
    q_points = range(
        -q_bz_ws_boundary + delta_q_for_sum / 2.0,
        q_bz_ws_boundary - delta_q_for_sum / 2.0,
        length = num_q_points
    )

    sum_integrand_complex = complex(0.0, 0.0)

    for q_val in q_points
        sum_integrand_complex += sigma_FM_integrand_1D(
            spec_params, e_k_bare, q_val,
            k_scalar_Ainv, hw_eval_ev, hw_phonon_ev,
            temp_k, efermi_ev, h_delta_ev
        )
    end

    integral_approximation = sum_integrand_complex * delta_q_for_sum
    prefactor = g0_1D_sq_ev2 / (2 * π)
    
    return prefactor * integral_approximation
end

function calc_sigma_FM_1D(
    spec_params::Spectral1D,
    e_k_bare::Function,
    k_array_Ainv::AbstractVector, hw_eval_array_ev::AbstractVector, g0_1D_sq_ev2::Float64, hw_phonon_ev::Float64,
    temp_k::Float64, efermi_ev::Float64, num_q_points::Int,
    h_delta_ev::Float64
)
    println("Calculating 1D Fan-Migdal self-enengy sigma(k,E)...")
    start_time = time()

    nk = length(k_array_Ainv)
    ne = length(hw_eval_array_ev)
    sigma_k_E_mat = zeros(Complex{Float64}, ne, nk)

    for (k_idx, k_val) in enumerate(k_array_Ainv)
        if (k_idx % max(1, nk ÷ 10) == 0) || k_idx == nk
            @printf(" Calculating for k-point %d/%d...\n", k_idx, nk)
        end
        # 并行计算能量点循环
        Threads.@threads for e_idx in 1:ne
            e_val = hw_eval_array_ev[e_idx]
            sigma_k_E_mat[e_idx, k_idx] = calc_sigma_k_e_FM_1D(
                spec_params, e_k_bare, k_val, e_val, g0_1D_sq_ev2,
                hw_phonon_ev, temp_k, efermi_ev, num_q_points, h_delta_ev
            )
        end
    end
    
    end_time = time()
    @printf("Calculation finished in %.2f seconds.\n\n", end_time - start_time)

    return sigma_k_E_mat
end

export Spectral1D, calc_sigma_FM_1D # 导出结构体和主函数

end