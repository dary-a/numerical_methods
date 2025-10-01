from math import sin, cos, sqrt


def trig_sin(x):
    un = x
    final = 0
    k = 0
    x_squared = x * x
    while abs(un) > epsilon_trig:
        final += un
        k += 2
        un = -un * x_squared / k / (k + 1)
    return final


def trig_cos(x):
    un = 1
    final = 0
    k = 0
    x_squared = x * x
    while abs(un) > epsilon_trig:
        final += un
        k += 2
        un = -un * x_squared / (k - 1) / k
    return final


def my_sqrt(x):
    ans = x
    prev_ans = 0
    # print(abs(ans - prev_ans))
    while abs(ans - prev_ans) > epsilon_sqrt:
        prev_ans = ans
        ans = 1 / 2 * (prev_ans + x / prev_ans)
        # print(abs(ans - prev_ans), abs(sqrt(x) - ans))
    return ans


epsilon_trig = 2.403 * 1e-7
epsilon_sqrt = 1.71 * 1e-7
left_border = 0.2
right_border = 0.3
step = 0.01
count = 0
stop = left_border

W = 15
header = (
    f"{'x':>{W}} | {'my_sin':>{W}} | {'math.sin':>{W}} | {'Δsin':>{W}} | "
    f"{'my_cos':>{W}} | {'math.cos':>{W}} | {'Δcos':>{W}} | "
    f"{'my_sqrt':>{W}} | {'math.sqrt':>{W}} | {'Δsqrt':>{W}} | "
    f"{'my_z':>{W}} | {'z':>{W}} | {'Δz':>{W}}"
)

print(header)
print("-" * len(header))

while stop <= right_border + 1e-12:
    my_sin_value, sin_value = trig_sin(3 * stop + 0.1), sin(3 * stop + 0.1)
    my_cos_value, cos_value = trig_cos(2 * stop + 0.3), cos(2 * stop + 0.3)

    my_sqrt_value, sqrt_value = my_sqrt(1 + stop * stop), sqrt(1 + stop * stop)
    my_z_value, z_value = my_sqrt_value * (my_sin_value + my_cos_value), sqrt_value * (sin_value + cos_value)
    # print(abs(my_sqrt_value - sqrt_value))
    print(
        f"{stop:{W}.3f} | {my_sin_value:{W}.8f} | {sin_value:{W}.8f} | {abs(my_sin_value - sin_value):{W}.8e} | "
        f"{my_cos_value:{W}.8f} | {cos_value:{W}.8f} | {abs(my_cos_value - cos_value):{W}.8e} | "
        f"{my_sqrt_value:{W}.8f} | {sqrt_value:{W}.8f} | {abs(my_sqrt_value - sqrt_value):{W}.8e} | "
        f"{my_z_value:{W}.8f} | {z_value:{W}.8f} | {abs(my_z_value - z_value):{W}.8e}"
    )
    count += 1
    stop = left_border + count * step