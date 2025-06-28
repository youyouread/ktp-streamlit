import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root, minimize_scalar, brentq
from PIL import Image

# =============================
# 折射率计算函数（Sellmeier）
# =============================
def calc_refractive_index(lam, direction):
    lam_um = lam / 1000.0
    if direction == 'x':
        n_sq = 3.0065 + 0.03901 / (lam_um ** 2 - 0.04251) - 0.01327 * lam_um ** 2
    elif direction == 'y':
        n_sq = 3.0333 + 0.04154 / (lam_um ** 2 - 0.04547) - 0.01408 * lam_um ** 2
    elif direction == 'z':
        n_sq = 3.3134 + 0.05694 / (lam_um ** 2 - 0.05658) - 0.01682 * lam_um ** 2
    else:
        return 0.0
    return np.sqrt(n_sq)

# =============================
# 效率计算所需：折射率平方函数
# =============================
def n_x_sq(lam):
    return 3.0065 + 0.03901 / (lam ** 2 - 0.04251) - 0.01327 * lam ** 2

def n_y_sq(lam):
    return 3.0333 + 0.04154 / (lam ** 2 - 0.04547) - 0.01408 * lam ** 2

def n_z_sq(lam):
    return 3.3134 + 0.05694 / (lam ** 2 - 0.05658) - 0.01682 * lam ** 2

# =============================
# e光折射率计算（用于相位匹配）
# =============================
def n_e(lam, phi):
    lam_um = lam / 1000.0
    n_x = n_x_sq(lam_um)
    n_y = n_y_sq(lam_um)
    return np.sqrt(n_x * np.sin(phi)**2 + n_y * np.cos(phi)**2)

# 相位匹配函数
def phase_matching_eq(phi, lam):
    lam1 = lam
    lam2 = lam1 / 2
    n1_o = calc_refractive_index(lam1, 'z')
    n1_e_val = n_e(lam1, phi)
    n2_e_val = n_e(lam2, phi)
    return 0.5 * (n1_o + n1_e_val) - n2_e_val

# =============================
# 效率计算函数
# =============================
def n1_e_deg(phi_deg):
    phi_rad = np.radians(phi_deg)
    return np.sqrt(nx1_sq * np.sin(phi_rad)**2 + ny1_sq * np.cos(phi_rad)**2)

def n2_e_deg(phi_deg):
    phi_rad = np.radians(phi_deg)
    return np.sqrt(nx2_sq * np.sin(phi_rad)**2 + ny2_sq * np.cos(phi_rad)**2)

def delta_k(phi_deg, L=10e-3):
    n1_e_val = n1_e_deg(phi_deg)
    n2_e_val = n2_e_deg(phi_deg)
    return (4 * np.pi / (lambda1 * 1e-6)) * (n2_e_val - 0.5 * (n1_o + n1_e_val))

def shg_efficiency(phi_deg, L=10e-3):
    dk = delta_k(phi_deg, L)
    dkL2 = dk * L / 2
    if dkL2 == 0:
        return 1.0
    return (np.sin(dkL2) / dkL2) ** 2

def find_optimal_phi():
    result = minimize_scalar(lambda phi: abs(delta_k(phi)), bounds=(0, 90), method='bounded')
    return result.x

# =============================
# Streamlit 主应用
# =============================
def main():
    import matplotlib
    matplotlib.rcParams['font.sans-serif'] = ['Microsoft YaHei']
    matplotlib.rcParams['axes.unicode_minus'] = False

    st.title("KTP晶体倍频二类匹配计算")

    st.sidebar.markdown("### 计算说明")
    st.sidebar.markdown("1. **折射率计算**：使用Sellmeier方程计算KTP晶体在三个方向(x, y, z)的折射率")
    st.sidebar.markdown("2. **二类相位匹配条件：**")
    st.sidebar.latex(r"\frac{1}{2} \left[ n_1^o + n_1^e(\theta, \varphi) \right] = n_2^e(\theta, \varphi)")
    st.sidebar.markdown("3. **角度求解**：计算相位匹配角φ随波长变化（θ固定为90°）")

    lam_input = st.sidebar.number_input("输入入射光波长 (nm)", 800.0, 1200.0, 1064.0, 1.0)


    # 折射率展示
    st.subheader("折射率计算结果")
    col1, col2 = st.columns(2)
    with col1:
        st.markdown(f"**基频光 (λ = {lam_input:.1f} nm)**")
        st.write(f"n_x = {calc_refractive_index(lam_input, 'x'):.4f}")
        st.write(f"n_y = {calc_refractive_index(lam_input, 'y'):.4f}")
        st.write(f"n_z = {calc_refractive_index(lam_input, 'z'):.4f}")
    with col2:
        st.markdown(f"**倍频光 (λ = {lam_input/2:.1f} nm)**")
        st.write(f"n_x = {calc_refractive_index(lam_input/2, 'x'):.4f}")
        st.write(f"n_y = {calc_refractive_index(lam_input/2, 'y'):.4f}")
        st.write(f"n_z = {calc_refractive_index(lam_input/2, 'z'):.4f}")

    # 折射率图
    st.subheader("折射率随波长变化曲线")
    wavelengths = np.linspace(500, 1200, 200)
    n_x = [calc_refractive_index(lam, 'x') for lam in wavelengths]
    n_y = [calc_refractive_index(lam, 'y') for lam in wavelengths]
    n_z = [calc_refractive_index(lam, 'z') for lam in wavelengths]

    fig1, ax1 = plt.subplots(figsize=(10, 6))
    ax1.plot(wavelengths, n_x, label='n_x', color='blue')
    ax1.plot(wavelengths, n_y, label='n_y', color='red')
    ax1.plot(wavelengths, n_z, label='n_z', color='darkgreen')
    ax1.axvline(lam_input, color='orange', linestyle='--', label=f'{lam_input}nm')
    ax1.axvline(lam_input / 2, color='purple', linestyle='--', label=f'{lam_input/2}nm')
    ax1.set_xlabel("波长 (nm)")
    ax1.set_ylabel("折射率")
    ax1.set_title("KTP晶体折射率 vs 波长")
    ax1.legend(loc='center left', bbox_to_anchor=(1.02, 0.5))
    ax1.grid(True)
    st.pyplot(fig1)

    # 相位匹配角图
    st.subheader("相位匹配角φ vs 波长 (1000–1100nm)")
    lam_range = np.linspace(1000, 1100, 300)
    def solve_phi_brentq(lam):
        try:
            return np.degrees(brentq(phase_matching_eq, np.radians(1), np.radians(89), args=(lam,)))
        except:
            return np.nan
    phi_vals = [solve_phi_brentq(l) for l in lam_range]

    fig2, ax2 = plt.subplots(figsize=(10, 6))
    ax2.plot(lam_range, phi_vals, 'b-')
    ax2.set_xlabel("波长 (nm)")
    ax2.set_ylabel("匹配角 φ (度)")
    ax2.set_ylim(0, 90)
    ax2.set_title("相位匹配角 vs 波长")
    ax2.grid(True)
    st.pyplot(fig2)

    # 效率计算
    st.subheader("倍频转换效率 vs 方位角 φ")
    global lambda1, lambda2, nx1_sq, ny1_sq, nz1_sq, nx2_sq, ny2_sq, nz2_sq, n1_o
    lambda1 = lam_input / 1000
    lambda2 = lambda1 / 2
    nx1_sq, ny1_sq, nz1_sq = n_x_sq(lambda1), n_y_sq(lambda1), n_z_sq(lambda1)
    nx2_sq, ny2_sq, nz2_sq = n_x_sq(lambda2), n_y_sq(lambda2), n_z_sq(lambda2)
    n1_o = np.sqrt(nz1_sq)

    phi_arr = np.linspace(0, 90, 500)
    eta_arr = [shg_efficiency(phi) for phi in phi_arr]
    opt_phi = find_optimal_phi()

    fig3, ax3 = plt.subplots(figsize=(10, 6))
    ax3.plot(phi_arr, eta_arr, 'b-', label="转换效率")
    ax3.axvline(opt_phi, color='r', linestyle='--', label=f"最佳角度: {opt_phi:.2f}°")
    ax3.set_xlabel("方位角φ (度)")
    ax3.set_ylabel("相对输出功率")
    ax3.set_ylim(0, 1.05)
    ax3.grid(True)
    ax3.legend()
    st.pyplot(fig3)

    st.sidebar.subheader("输入φ角计算效率")
    phi_sel = st.sidebar.number_input("输入φ角 (度)", 0.0, 90.0, 24.93, 1.0)


    # 放晶体图片在最后
    try:
        image = Image.open("ktp.png")
        st.sidebar.image(image, caption="KTP晶体示意图", use_container_width=True)
    except:
        st.sidebar.warning("未找到晶体图像")
    st.write(f"**当前选择角度：{phi_sel:.2f}°**")
    st.write(f"- 转换效率 η = {shg_efficiency(phi_sel):.6f}")
    st.write(f"- Δk = {delta_k(phi_sel):.3e} rad/m")
    st.write(f"- n1⁽ᵉ⁾ = {n1_e_deg(phi_sel):.6f}")
    st.write(f"- n2⁽ᵉ⁾ = {n2_e_deg(phi_sel):.6f}")
    st.write(f"- n1⁽ᵒ⁾ = {n1_o:.6f}")



if __name__ == "__main__":
    main()
