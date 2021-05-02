
from derivatives import Derivative, print_dict, print_dict_as_code_cpp_format


if __name__ == "__main__":


    du2_normal, dv2_normal = Derivative(small_angles=False).compute_derivative_2nd()
    du2_small, dv2_small = Derivative(small_angles=True).compute_derivative_2nd()

    print("")
    print("/***************************************************************/")
    print("/********************** Derivative U ***************************/")
    print("/***************************************************************/")
    print(print_dict_as_code_cpp_format(du2_normal, du2_small))
    print("")

    print("/***************************************************************/")
    print("/********************** Derivative V ***************************/")
    print("/***************************************************************/")
    print(print_dict_as_code_cpp_format(dv2_normal, dv2_small))
    print("")
