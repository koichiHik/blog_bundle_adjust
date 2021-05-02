
from derivatives import Derivative, print_dict, print_dict_as_code_cpp_format


if __name__ == "__main__":


    du_normal, dv_normal = Derivative(small_angles=False).compute_derivative_1st()
    du_small, dv_small = Derivative(small_angles=True).compute_derivative_1st()

    print("")
    print("/***************************************************************/")
    print("/********************** Derivative U ***************************/")
    print("/***************************************************************/")
    print(print_dict_as_code_cpp_format(du_normal, du_small))
    print("")

    print("/***************************************************************/")
    print("/********************** Derivative V ***************************/")
    print("/***************************************************************/")
    print(print_dict_as_code_cpp_format(dv_normal, dv_small))
    print("")
