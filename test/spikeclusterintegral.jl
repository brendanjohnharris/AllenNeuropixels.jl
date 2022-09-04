using SymPy

t, t₁, t₂, σ = symbols("t, t₁, t₂, σ", real=true)


f = simplify((exp((-1/2)*((t - t₁)/σ)^2)*exp((-1/2)*((t - t₂)/σ)^2))/((σ*sqrt(PI))*(σ*sqrt(2PI))))
I = integrate(f, (t, -oo, oo)) |> simplify
SymPy.latex(I)

# \begin{cases} \frac{0.5 \sqrt{2} e^{\frac{- 0.25 t_1^{2} + 0.5 t_1 t_2 - 0.25 t_2^{2}}{s^{2}}}}{\sqrt{\pi} s} & \text{for}\: \left(\left|{\arg{\left(s \right)}}\right| \leq \frac{\pi}{4} \wedge \left|{2 \arg{\left(t_1 \right)} - 4 \arg{\left(s \right)} + 2 \arg{\left(\frac{t_1 - t_2}{t_1} \right)} + 2 \pi}\right| < \pi \wedge \left|{2 \arg{\left(t_2 \right)} - 4 \arg{\left(s \right)} + 2 \arg{\left(\frac{- t_1 + t_2}{t_2} \right)} + 2 \pi}\right| < \pi\right) \vee \left|{\arg{\left(s \right)}}\right| < \frac{\pi}{4} \\\frac{\sqrt{2} \int\limits_{-\infty}^{\infty} e^{- \frac{0.5 \left(\left(t - t_1\right)^{2} + \left(t - t_2\right)^{2}\right)}{s^{2}}}\, dt}{2 \pi s^{2}} & \text{otherwise} \end{cases}
