
import sympy
from sympy import ccode
from sympy.printing import print_ccode
from collections import OrderedDict


def print_dict(derive_dict):

    for k, v in derive_dict.items():
        print(k)
        print("")
        print("return " + str(v) + ";")
        print("")

def print_dict_as_code_cpp_format(derive_dict_normal, derive_dict_small):

    lines = ''
    for k, v in derive_dict_normal.items():
        lines = lines \
            + "\n\ninline double " + str(k) \
            + "(const Eigen::Vector3d &T, const Eigen::Vector3d &rot," \
            + "const Eigen::Matrix3d &K, const Eigen::Vector3d &p, const Eigen::Vector3d &x) {" \
            + "double theta = std::pow(std::pow(rot(0),2) + std::pow(rot(1),2) + std::pow(rot(2),2), 0.5);" \
            + "if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {" \
            + "return " + str(ccode(v, standard='C89')) + ";" \
            + "} else {" \
            + "return " + str(ccode(derive_dict_small[k], standard='C89')) + ";" \
            + "}" \
            + "}"

    return lines

class Derivative:

    def __init__(self, small_angles=False):

        # 1. 3D Point Coordinate.
        self.X = sympy.Symbol('p(0)')
        self.Y = sympy.Symbol('p(1)')
        self.Z = sympy.Symbol('p(2)')

        # 2. Internal Camera Parameters.
        self.fx = sympy.Symbol('K(0,0)')
        self.fy = sympy.Symbol('K(1,1)')
        self.cx = sympy.Symbol('K(0,2)')
        self.cy = sympy.Symbol('K(1,2)')
        self.s = sympy.Symbol('K(0,1)')
        self.f0 = sympy.Symbol('K(2,2)')
        self.K = sympy.Matrix([[self.fx, self.s, self.cx], [0, self.fy, self.cy], [0, 0, self.f0]])

        # 3. External Camera Parameters.
        # Translation.
        self.tx = sympy.Symbol('T(0)')
        self.ty = sympy.Symbol('T(1)')
        self.tz = sympy.Symbol('T(2)')
        self.T = sympy.Matrix([self.tx, self.ty, self.tz])

        # Rotation.
        self.vx = sympy.Symbol('rot(0)')
        self.vy = sympy.Symbol('rot(1)')
        self.vz = sympy.Symbol('rot(2)')
        self.theta = sympy.sqrt(self.vx**2 + self.vy**2 + self.vz**2)


        self.I = sympy.Matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        self.W = sympy.Matrix([[0, -self.vz, self.vy], [self.vz, 0, -self.vx], [-self.vy, self.vx, 0]])

        if (not small_angles):
            self.abbre = self.abbre_normal_angle
            self.R = self.I + (sympy.sin(self.theta) / self.theta) * self.W + \
                ((1 - sympy.cos(self.theta)) / self.theta**2) * self.W * self.W
            self.P = self.K * (self.R.row_join(self.T))
        else:
            self.abbre = self.abbre_small_angle
            self.R = self.I + self.W
            self.P = self.K * (self.R.row_join(self.T))

        self.x = self.P * sympy.Matrix([self.X, self.Y, self.Z, 1])
        self.u = self.x[0] / self.x[2]
        self.v = self.x[1] / self.x[2]


    def compute_derivative_1st(self, abbreviate=True):

        if (not abbreviate):
            def empty_abbreviate(expr):
                return expr
            abbre = empty_abbreviate
        else:
            abbre = self.abbre

        # Derivative regarding u.
        du = OrderedDict()
        du["du_dX"] = abbre(sympy.diff(self.u, self.X))
        du["du_dY"] = abbre(sympy.diff(self.u, self.Y))
        du["du_dZ"] = abbre(sympy.diff(self.u, self.Z))
        du["du_dfx"] = abbre(sympy.diff(self.u, self.fx))
        du["du_dfy"] = abbre(sympy.diff(self.u, self.fy))
        du["du_dcx"] = abbre(sympy.diff(self.u, self.cx))
        du["du_dcy"] = abbre(sympy.diff(self.u, self.cy))
        du["du_ds"] = abbre(sympy.diff(self.u, self.s))
        du["du_dtx"] = abbre(sympy.diff(self.u, self.tx))
        du["du_dty"] = abbre(sympy.diff(self.u, self.ty))
        du["du_dtz"] = abbre(sympy.diff(self.u, self.tz))
        du["du_dvx"] = abbre(sympy.diff(self.u, self.vx))
        du["du_dvy"] = abbre(sympy.diff(self.u, self.vy))
        du["du_dvz"] = abbre(sympy.diff(self.u, self.vz))

        # Derivative regarding v.
        dv = OrderedDict()
        dv["dv_dX"] = abbre(sympy.diff(self.v, self.X))
        dv["dv_dY"] = abbre(sympy.diff(self.v, self.Y))
        dv["dv_dZ"] = abbre(sympy.diff(self.v, self.Z))
        dv["dv_dfx"] = abbre(sympy.diff(self.v, self.fx))
        dv["dv_dfy"] = abbre(sympy.diff(self.v, self.fy))
        dv["dv_dcx"] = abbre(sympy.diff(self.v, self.cx))
        dv["dv_dcy"] = abbre(sympy.diff(self.v, self.cy))
        dv["dv_ds"] = abbre(sympy.diff(self.v, self.s))
        dv["dv_dtx"] = abbre(sympy.diff(self.v, self.tx))
        dv["dv_dty"] = abbre(sympy.diff(self.v, self.ty))
        dv["dv_dtz"] = abbre(sympy.diff(self.v, self.tz))
        dv["dv_dvx"] = abbre(sympy.diff(self.v, self.vx))
        dv["dv_dvy"] = abbre(sympy.diff(self.v, self.vy))
        dv["dv_dvz"] = abbre(sympy.diff(self.v, self.vz))

        return du, dv


    def compute_derivative_2nd(self):

        du, dv = self.compute_derivative_1st(False)

        # Derivative regarding u.
        du2 = OrderedDict()
        du2["du2_dXdX"] = self.abbre(sympy.diff(du["du_dX"], self.X))
        du2["du2_dXdY"] = self.abbre(sympy.diff(du["du_dX"], self.Y))
        du2["du2_dXdZ"] = self.abbre(sympy.diff(du["du_dX"], self.Z))
        du2["du2_dXdfx"] = self.abbre(sympy.diff(du["du_dX"], self.fx))
        du2["du2_dXdfy"] = self.abbre(sympy.diff(du["du_dX"], self.fy))
        du2["du2_dXdcx"] = self.abbre(sympy.diff(du["du_dX"], self.cx))
        du2["du2_dXdcy"] = self.abbre(sympy.diff(du["du_dX"], self.cy))
        du2["du2_dXds"] = self.abbre(sympy.diff(du["du_dX"], self.s))
        du2["du2_dXdtx"] = self.abbre(sympy.diff(du["du_dX"], self.tx))
        du2["du2_dXdty"] = self.abbre(sympy.diff(du["du_dX"], self.ty))
        du2["du2_dXdtz"] = self.abbre(sympy.diff(du["du_dX"], self.tz))
        du2["du2_dXdvx"] = self.abbre(sympy.diff(du["du_dX"], self.vx))
        du2["du2_dXdvy"] = self.abbre(sympy.diff(du["du_dX"], self.vy))
        du2["du2_dXdvz"] = self.abbre(sympy.diff(du["du_dX"], self.vz))

        du2["du2_dYdX"] = self.abbre(sympy.diff(du["du_dY"], self.X))
        du2["du2_dYdY"] = self.abbre(sympy.diff(du["du_dY"], self.Y))
        du2["du2_dYdZ"] = self.abbre(sympy.diff(du["du_dY"], self.Z))
        du2["du2_dYdfx"] = self.abbre(sympy.diff(du["du_dY"], self.fx))
        du2["du2_dYdfy"] = self.abbre(sympy.diff(du["du_dY"], self.fy))
        du2["du2_dYdcx"] = self.abbre(sympy.diff(du["du_dY"], self.cx))
        du2["du2_dYdcy"] = self.abbre(sympy.diff(du["du_dY"], self.cy))
        du2["du2_dYds"] = self.abbre(sympy.diff(du["du_dY"], self.s))
        du2["du2_dYdtx"] = self.abbre(sympy.diff(du["du_dY"], self.tx))
        du2["du2_dYdty"] = self.abbre(sympy.diff(du["du_dY"], self.ty))
        du2["du2_dYdtz"] = self.abbre(sympy.diff(du["du_dY"], self.tz))
        du2["du2_dYdvx"] = self.abbre(sympy.diff(du["du_dY"], self.vx))
        du2["du2_dYdvy"] = self.abbre(sympy.diff(du["du_dY"], self.vy))
        du2["du2_dYdvz"] = self.abbre(sympy.diff(du["du_dY"], self.vz))

        du2["du2_dZdX"] = self.abbre(sympy.diff(du["du_dZ"], self.X))
        du2["du2_dZdY"] = self.abbre(sympy.diff(du["du_dZ"], self.Y))
        du2["du2_dZdZ"] = self.abbre(sympy.diff(du["du_dZ"], self.Z))
        du2["du2_dZdfx"] = self.abbre(sympy.diff(du["du_dZ"], self.fx))
        du2["du2_dZdfy"] = self.abbre(sympy.diff(du["du_dZ"], self.fy))
        du2["du2_dZdcx"] = self.abbre(sympy.diff(du["du_dZ"], self.cx))
        du2["du2_dZdcy"] = self.abbre(sympy.diff(du["du_dZ"], self.cy))
        du2["du2_dZds"] = self.abbre(sympy.diff(du["du_dZ"], self.s))
        du2["du2_dZdtx"] = self.abbre(sympy.diff(du["du_dZ"], self.tx))
        du2["du2_dZdty"] = self.abbre(sympy.diff(du["du_dZ"], self.ty))
        du2["du2_dZdtz"] = self.abbre(sympy.diff(du["du_dZ"], self.tz))
        du2["du2_dZdvx"] = self.abbre(sympy.diff(du["du_dZ"], self.vx))
        du2["du2_dZdvy"] = self.abbre(sympy.diff(du["du_dZ"], self.vy))
        du2["du2_dZdvz"] = self.abbre(sympy.diff(du["du_dZ"], self.vz))

        du2["du2_dfxdX"] = self.abbre(sympy.diff(du["du_dfx"], self.X))
        du2["du2_dfxdY"] = self.abbre(sympy.diff(du["du_dfx"], self.Y))
        du2["du2_dfxdZ"] = self.abbre(sympy.diff(du["du_dfx"], self.Z))
        du2["du2_dfxdfx"] = self.abbre(sympy.diff(du["du_dfx"], self.fx))
        du2["du2_dfxdfy"] = self.abbre(sympy.diff(du["du_dfx"], self.fy))
        du2["du2_dfxdcx"] = self.abbre(sympy.diff(du["du_dfx"], self.cx))
        du2["du2_dfxdcy"] = self.abbre(sympy.diff(du["du_dfx"], self.cy))
        du2["du2_dfxds"] = self.abbre(sympy.diff(du["du_dfx"], self.s))
        du2["du2_dfxdtx"] = self.abbre(sympy.diff(du["du_dfx"], self.tx))
        du2["du2_dfxdty"] = self.abbre(sympy.diff(du["du_dfx"], self.ty))
        du2["du2_dfxdtz"] = self.abbre(sympy.diff(du["du_dfx"], self.tz))
        du2["du2_dfxdvx"] = self.abbre(sympy.diff(du["du_dfx"], self.vx))
        du2["du2_dfxdvy"] = self.abbre(sympy.diff(du["du_dfx"], self.vy))
        du2["du2_dfxdvz"] = self.abbre(sympy.diff(du["du_dfx"], self.vz))

        du2["du2_dfydX"] = self.abbre(sympy.diff(du["du_dfy"], self.X))
        du2["du2_dfydY"] = self.abbre(sympy.diff(du["du_dfy"], self.Y))
        du2["du2_dfydZ"] = self.abbre(sympy.diff(du["du_dfy"], self.Z))
        du2["du2_dfydfx"] = self.abbre(sympy.diff(du["du_dfy"], self.fx))
        du2["du2_dfydfy"] = self.abbre(sympy.diff(du["du_dfy"], self.fy))
        du2["du2_dfydcx"] = self.abbre(sympy.diff(du["du_dfy"], self.cx))
        du2["du2_dfydcy"] = self.abbre(sympy.diff(du["du_dfy"], self.cy))
        du2["du2_dfyds"] = self.abbre(sympy.diff(du["du_dfy"], self.s))
        du2["du2_dfydtx"] = self.abbre(sympy.diff(du["du_dfy"], self.tx))
        du2["du2_dfydty"] = self.abbre(sympy.diff(du["du_dfy"], self.ty))
        du2["du2_dfydtz"] = self.abbre(sympy.diff(du["du_dfy"], self.tz))
        du2["du2_dfydvx"] = self.abbre(sympy.diff(du["du_dfy"], self.vx))
        du2["du2_dfydvy"] = self.abbre(sympy.diff(du["du_dfy"], self.vy))
        du2["du2_dfydvz"] = self.abbre(sympy.diff(du["du_dfy"], self.vz))

        du2["du2_dcxdX"] = self.abbre(sympy.diff(du["du_dcx"], self.X))
        du2["du2_dcxdY"] = self.abbre(sympy.diff(du["du_dcx"], self.Y))
        du2["du2_dcxdZ"] = self.abbre(sympy.diff(du["du_dcx"], self.Z))
        du2["du2_dcxdfx"] = self.abbre(sympy.diff(du["du_dcx"], self.fx))
        du2["du2_dcxdfy"] = self.abbre(sympy.diff(du["du_dcx"], self.fy))
        du2["du2_dcxdcx"] = self.abbre(sympy.diff(du["du_dcx"], self.cx))
        du2["du2_dcxdcy"] = self.abbre(sympy.diff(du["du_dcx"], self.cy))
        du2["du2_dcxds"] = self.abbre(sympy.diff(du["du_dcx"], self.s))
        du2["du2_dcxdtx"] = self.abbre(sympy.diff(du["du_dcx"], self.tx))
        du2["du2_dcxdty"] = self.abbre(sympy.diff(du["du_dcx"], self.ty))
        du2["du2_dcxdtz"] = self.abbre(sympy.diff(du["du_dcx"], self.tz))
        du2["du2_dcxdvx"] = self.abbre(sympy.diff(du["du_dcx"], self.vx))
        du2["du2_dcxdvy"] = self.abbre(sympy.diff(du["du_dcx"], self.vy))
        du2["du2_dcxdvz"] = self.abbre(sympy.diff(du["du_dcx"], self.vz))

        du2["du2_dcydX"] = self.abbre(sympy.diff(du["du_dcy"], self.X))
        du2["du2_dcydY"] = self.abbre(sympy.diff(du["du_dcy"], self.Y))
        du2["du2_dcydZ"] = self.abbre(sympy.diff(du["du_dcy"], self.Z))
        du2["du2_dcydfx"] = self.abbre(sympy.diff(du["du_dcy"], self.fx))
        du2["du2_dcydfy"] = self.abbre(sympy.diff(du["du_dcy"], self.fy))
        du2["du2_dcydcx"] = self.abbre(sympy.diff(du["du_dcy"], self.cx))
        du2["du2_dcydcy"] = self.abbre(sympy.diff(du["du_dcy"], self.cy))
        du2["du2_dcyds"] = self.abbre(sympy.diff(du["du_dcy"], self.s))
        du2["du2_dcydtx"] = self.abbre(sympy.diff(du["du_dcy"], self.tx))
        du2["du2_dcydty"] = self.abbre(sympy.diff(du["du_dcy"], self.ty))
        du2["du2_dcydtz"] = self.abbre(sympy.diff(du["du_dcy"], self.tz))
        du2["du2_dcydvx"] = self.abbre(sympy.diff(du["du_dcy"], self.vx))
        du2["du2_dcydvy"] = self.abbre(sympy.diff(du["du_dcy"], self.vy))
        du2["du2_dcydvz"] = self.abbre(sympy.diff(du["du_dcy"], self.vz))

        du2["du2_dsdX"] = self.abbre(sympy.diff(du["du_ds"], self.X))
        du2["du2_dsdY"] = self.abbre(sympy.diff(du["du_ds"], self.Y))
        du2["du2_dsdZ"] = self.abbre(sympy.diff(du["du_ds"], self.Z))
        du2["du2_dsdfx"] = self.abbre(sympy.diff(du["du_ds"], self.fx))
        du2["du2_dsdfy"] = self.abbre(sympy.diff(du["du_ds"], self.fy))
        du2["du2_dsdcx"] = self.abbre(sympy.diff(du["du_ds"], self.cx))
        du2["du2_dsdcy"] = self.abbre(sympy.diff(du["du_ds"], self.cy))
        du2["du2_dsds"] = self.abbre(sympy.diff(du["du_ds"], self.s))
        du2["du2_dsdtx"] = self.abbre(sympy.diff(du["du_ds"], self.tx))
        du2["du2_dsdty"] = self.abbre(sympy.diff(du["du_ds"], self.ty))
        du2["du2_dsdtz"] = self.abbre(sympy.diff(du["du_ds"], self.tz))
        du2["du2_dsdvx"] = self.abbre(sympy.diff(du["du_ds"], self.vx))
        du2["du2_dsdvy"] = self.abbre(sympy.diff(du["du_ds"], self.vy))
        du2["du2_dsdvz"] = self.abbre(sympy.diff(du["du_ds"], self.vz))

        du2["du2_dtxdX"] = self.abbre(sympy.diff(du["du_dtx"], self.X))
        du2["du2_dtxdY"] = self.abbre(sympy.diff(du["du_dtx"], self.Y))
        du2["du2_dtxdZ"] = self.abbre(sympy.diff(du["du_dtx"], self.Z))
        du2["du2_dtxdfx"] = self.abbre(sympy.diff(du["du_dtx"], self.fx))
        du2["du2_dtxdfy"] = self.abbre(sympy.diff(du["du_dtx"], self.fy))
        du2["du2_dtxdcx"] = self.abbre(sympy.diff(du["du_dtx"], self.cx))
        du2["du2_dtxdcy"] = self.abbre(sympy.diff(du["du_dtx"], self.cy))
        du2["du2_dtxds"] = self.abbre(sympy.diff(du["du_dtx"], self.s))
        du2["du2_dtxdtx"] = self.abbre(sympy.diff(du["du_dtx"], self.tx))
        du2["du2_dtxdty"] = self.abbre(sympy.diff(du["du_dtx"], self.ty))
        du2["du2_dtxdtz"] = self.abbre(sympy.diff(du["du_dtx"], self.tz))
        du2["du2_dtxdvx"] = self.abbre(sympy.diff(du["du_dtx"], self.vx))
        du2["du2_dtxdvy"] = self.abbre(sympy.diff(du["du_dtx"], self.vy))
        du2["du2_dtxdvz"] = self.abbre(sympy.diff(du["du_dtx"], self.vz))

        du2["du2_dtydX"] = self.abbre(sympy.diff(du["du_dty"], self.X))
        du2["du2_dtydY"] = self.abbre(sympy.diff(du["du_dty"], self.Y))
        du2["du2_dtydZ"] = self.abbre(sympy.diff(du["du_dty"], self.Z))
        du2["du2_dtydfx"] = self.abbre(sympy.diff(du["du_dty"], self.fx))
        du2["du2_dtydfy"] = self.abbre(sympy.diff(du["du_dty"], self.fy))
        du2["du2_dtydcx"] = self.abbre(sympy.diff(du["du_dty"], self.cx))
        du2["du2_dtydcy"] = self.abbre(sympy.diff(du["du_dty"], self.cy))
        du2["du2_dtyds"] = self.abbre(sympy.diff(du["du_dty"], self.s))
        du2["du2_dtydtx"] = self.abbre(sympy.diff(du["du_dty"], self.tx))
        du2["du2_dtydty"] = self.abbre(sympy.diff(du["du_dty"], self.ty))
        du2["du2_dtydtz"] = self.abbre(sympy.diff(du["du_dty"], self.tz))
        du2["du2_dtydvx"] = self.abbre(sympy.diff(du["du_dty"], self.vx))
        du2["du2_dtydvy"] = self.abbre(sympy.diff(du["du_dty"], self.vy))
        du2["du2_dtydvz"] = self.abbre(sympy.diff(du["du_dty"], self.vz))

        du2["du2_dtzdX"] = self.abbre(sympy.diff(du["du_dtz"], self.X))
        du2["du2_dtzdY"] = self.abbre(sympy.diff(du["du_dtz"], self.Y))
        du2["du2_dtzdZ"] = self.abbre(sympy.diff(du["du_dtz"], self.Z))
        du2["du2_dtzdfx"] = self.abbre(sympy.diff(du["du_dtz"], self.fx))
        du2["du2_dtzdfy"] = self.abbre(sympy.diff(du["du_dtz"], self.fy))
        du2["du2_dtzdcx"] = self.abbre(sympy.diff(du["du_dtz"], self.cx))
        du2["du2_dtzdcy"] = self.abbre(sympy.diff(du["du_dtz"], self.cy))
        du2["du2_dtzds"] = self.abbre(sympy.diff(du["du_dtz"], self.s))
        du2["du2_dtzdtx"] = self.abbre(sympy.diff(du["du_dtz"], self.tx))
        du2["du2_dtzdty"] = self.abbre(sympy.diff(du["du_dtz"], self.ty))
        du2["du2_dtzdtz"] = self.abbre(sympy.diff(du["du_dtz"], self.tz))
        du2["du2_dtzdvx"] = self.abbre(sympy.diff(du["du_dtz"], self.vx))
        du2["du2_dtzdvy"] = self.abbre(sympy.diff(du["du_dtz"], self.vy))
        du2["du2_dtzdvz"] = self.abbre(sympy.diff(du["du_dtz"], self.vz))

        du2["du2_dvxdX"] = self.abbre(sympy.diff(du["du_dvx"], self.X))
        du2["du2_dvxdY"] = self.abbre(sympy.diff(du["du_dvx"], self.Y))
        du2["du2_dvxdZ"] = self.abbre(sympy.diff(du["du_dvx"], self.Z))
        du2["du2_dvxdfx"] = self.abbre(sympy.diff(du["du_dvx"], self.fx))
        du2["du2_dvxdfy"] = self.abbre(sympy.diff(du["du_dvx"], self.fy))
        du2["du2_dvxdcx"] = self.abbre(sympy.diff(du["du_dvx"], self.cx))
        du2["du2_dvxdcy"] = self.abbre(sympy.diff(du["du_dvx"], self.cy))
        du2["du2_dvxds"] = self.abbre(sympy.diff(du["du_dvx"], self.s))
        du2["du2_dvxdtx"] = self.abbre(sympy.diff(du["du_dvx"], self.tx))
        du2["du2_dvxdty"] = self.abbre(sympy.diff(du["du_dvx"], self.ty))
        du2["du2_dvxdtz"] = self.abbre(sympy.diff(du["du_dvx"], self.tz))
        du2["du2_dvxdvx"] = self.abbre(sympy.diff(du["du_dvx"], self.vx))
        du2["du2_dvxdvy"] = self.abbre(sympy.diff(du["du_dvx"], self.vy))
        du2["du2_dvxdvz"] = self.abbre(sympy.diff(du["du_dvx"], self.vz))

        du2["du2_dvydX"] = self.abbre(sympy.diff(du["du_dvy"], self.X))
        du2["du2_dvydY"] = self.abbre(sympy.diff(du["du_dvy"], self.Y))
        du2["du2_dvydZ"] = self.abbre(sympy.diff(du["du_dvy"], self.Z))
        du2["du2_dvydfx"] = self.abbre(sympy.diff(du["du_dvy"], self.fx))
        du2["du2_dvydfy"] = self.abbre(sympy.diff(du["du_dvy"], self.fy))
        du2["du2_dvydcx"] = self.abbre(sympy.diff(du["du_dvy"], self.cx))
        du2["du2_dvydcy"] = self.abbre(sympy.diff(du["du_dvy"], self.cy))
        du2["du2_dvyds"] = self.abbre(sympy.diff(du["du_dvy"], self.s))
        du2["du2_dvydtx"] = self.abbre(sympy.diff(du["du_dvy"], self.tx))
        du2["du2_dvydty"] = self.abbre(sympy.diff(du["du_dvy"], self.ty))
        du2["du2_dvydtz"] = self.abbre(sympy.diff(du["du_dvy"], self.tz))
        du2["du2_dvydvx"] = self.abbre(sympy.diff(du["du_dvy"], self.vx))
        du2["du2_dvydvy"] = self.abbre(sympy.diff(du["du_dvy"], self.vy))
        du2["du2_dvydvz"] = self.abbre(sympy.diff(du["du_dvy"], self.vz))

        du2["du2_dvzdX"] = self.abbre(sympy.diff(du["du_dvz"], self.X))
        du2["du2_dvzdY"] = self.abbre(sympy.diff(du["du_dvz"], self.Y))
        du2["du2_dvzdZ"] = self.abbre(sympy.diff(du["du_dvz"], self.Z))
        du2["du2_dvzdfx"] = self.abbre(sympy.diff(du["du_dvz"], self.fx))
        du2["du2_dvzdfy"] = self.abbre(sympy.diff(du["du_dvz"], self.fy))
        du2["du2_dvzdcx"] = self.abbre(sympy.diff(du["du_dvz"], self.cx))
        du2["du2_dvzdcy"] = self.abbre(sympy.diff(du["du_dvz"], self.cy))
        du2["du2_dvzds"] = self.abbre(sympy.diff(du["du_dvz"], self.s))
        du2["du2_dvzdtx"] = self.abbre(sympy.diff(du["du_dvz"], self.tx))
        du2["du2_dvzdty"] = self.abbre(sympy.diff(du["du_dvz"], self.ty))
        du2["du2_dvzdtz"] = self.abbre(sympy.diff(du["du_dvz"], self.tz))
        du2["du2_dvzdvx"] = self.abbre(sympy.diff(du["du_dvz"], self.vx))
        du2["du2_dvzdvy"] = self.abbre(sympy.diff(du["du_dvz"], self.vy))
        du2["du2_dvzdvz"] = self.abbre(sympy.diff(du["du_dvz"], self.vz))


        # Derivative regarding v.
        dv2 = OrderedDict()
        dv2["dv2_dXdX"] = self.abbre(sympy.diff(dv["dv_dX"], self.X))
        dv2["dv2_dXdY"] = self.abbre(sympy.diff(dv["dv_dX"], self.Y))
        dv2["dv2_dXdZ"] = self.abbre(sympy.diff(dv["dv_dX"], self.Z))
        dv2["dv2_dXdfx"] = self.abbre(sympy.diff(dv["dv_dX"], self.fx))
        dv2["dv2_dXdfy"] = self.abbre(sympy.diff(dv["dv_dX"], self.fy))
        dv2["dv2_dXdcx"] = self.abbre(sympy.diff(dv["dv_dX"], self.cx))
        dv2["dv2_dXdcy"] = self.abbre(sympy.diff(dv["dv_dX"], self.cy))
        dv2["dv2_dXds"] = self.abbre(sympy.diff(dv["dv_dX"], self.s))
        dv2["dv2_dXdtx"] = self.abbre(sympy.diff(dv["dv_dX"], self.tx))
        dv2["dv2_dXdty"] = self.abbre(sympy.diff(dv["dv_dX"], self.ty))
        dv2["dv2_dXdtz"] = self.abbre(sympy.diff(dv["dv_dX"], self.tz))
        dv2["dv2_dXdvx"] = self.abbre(sympy.diff(dv["dv_dX"], self.vx))
        dv2["dv2_dXdvy"] = self.abbre(sympy.diff(dv["dv_dX"], self.vy))
        dv2["dv2_dXdvz"] = self.abbre(sympy.diff(dv["dv_dX"], self.vz))

        dv2["dv2_dYdX"] = self.abbre(sympy.diff(dv["dv_dY"], self.X))
        dv2["dv2_dYdY"] = self.abbre(sympy.diff(dv["dv_dY"], self.Y))
        dv2["dv2_dYdZ"] = self.abbre(sympy.diff(dv["dv_dY"], self.Z))
        dv2["dv2_dYdfx"] = self.abbre(sympy.diff(dv["dv_dY"], self.fx))
        dv2["dv2_dYdfy"] = self.abbre(sympy.diff(dv["dv_dY"], self.fy))
        dv2["dv2_dYdcx"] = self.abbre(sympy.diff(dv["dv_dY"], self.cx))
        dv2["dv2_dYdcy"] = self.abbre(sympy.diff(dv["dv_dY"], self.cy))
        dv2["dv2_dYds"] = self.abbre(sympy.diff(dv["dv_dY"], self.s))
        dv2["dv2_dYdtx"] = self.abbre(sympy.diff(dv["dv_dY"], self.tx))
        dv2["dv2_dYdty"] = self.abbre(sympy.diff(dv["dv_dY"], self.ty))
        dv2["dv2_dYdtz"] = self.abbre(sympy.diff(dv["dv_dY"], self.tz))
        dv2["dv2_dYdvx"] = self.abbre(sympy.diff(dv["dv_dY"], self.vx))
        dv2["dv2_dYdvy"] = self.abbre(sympy.diff(dv["dv_dY"], self.vy))
        dv2["dv2_dYdvz"] = self.abbre(sympy.diff(dv["dv_dY"], self.vz))

        dv2["dv2_dZdX"] = self.abbre(sympy.diff(dv["dv_dZ"], self.X))
        dv2["dv2_dZdY"] = self.abbre(sympy.diff(dv["dv_dZ"], self.Y))
        dv2["dv2_dZdZ"] = self.abbre(sympy.diff(dv["dv_dZ"], self.Z))
        dv2["dv2_dZdfx"] = self.abbre(sympy.diff(dv["dv_dZ"], self.fx))
        dv2["dv2_dZdfy"] = self.abbre(sympy.diff(dv["dv_dZ"], self.fy))
        dv2["dv2_dZdcx"] = self.abbre(sympy.diff(dv["dv_dZ"], self.cx))
        dv2["dv2_dZdcy"] = self.abbre(sympy.diff(dv["dv_dZ"], self.cy))
        dv2["dv2_dZds"] = self.abbre(sympy.diff(dv["dv_dZ"], self.s))
        dv2["dv2_dZdtx"] = self.abbre(sympy.diff(dv["dv_dZ"], self.tx))
        dv2["dv2_dZdty"] = self.abbre(sympy.diff(dv["dv_dZ"], self.ty))
        dv2["dv2_dZdtz"] = self.abbre(sympy.diff(dv["dv_dZ"], self.tz))
        dv2["dv2_dZdvx"] = self.abbre(sympy.diff(dv["dv_dZ"], self.vx))
        dv2["dv2_dZdvy"] = self.abbre(sympy.diff(dv["dv_dZ"], self.vy))
        dv2["dv2_dZdvz"] = self.abbre(sympy.diff(dv["dv_dZ"], self.vz))

        dv2["dv2_dfxdX"] = self.abbre(sympy.diff(dv["dv_dfx"], self.X))
        dv2["dv2_dfxdY"] = self.abbre(sympy.diff(dv["dv_dfx"], self.Y))
        dv2["dv2_dfxdZ"] = self.abbre(sympy.diff(dv["dv_dfx"], self.Z))
        dv2["dv2_dfxdfx"] = self.abbre(sympy.diff(dv["dv_dfx"], self.fx))
        dv2["dv2_dfxdfy"] = self.abbre(sympy.diff(dv["dv_dfx"], self.fy))
        dv2["dv2_dfxdcx"] = self.abbre(sympy.diff(dv["dv_dfx"], self.cx))
        dv2["dv2_dfxdcy"] = self.abbre(sympy.diff(dv["dv_dfx"], self.cy))
        dv2["dv2_dfxds"] = self.abbre(sympy.diff(dv["dv_dfx"], self.s))
        dv2["dv2_dfxdtx"] = self.abbre(sympy.diff(dv["dv_dfx"], self.tx))
        dv2["dv2_dfxdty"] = self.abbre(sympy.diff(dv["dv_dfx"], self.ty))
        dv2["dv2_dfxdtz"] = self.abbre(sympy.diff(dv["dv_dfx"], self.tz))
        dv2["dv2_dfxdvx"] = self.abbre(sympy.diff(dv["dv_dfx"], self.vx))
        dv2["dv2_dfxdvy"] = self.abbre(sympy.diff(dv["dv_dfx"], self.vy))
        dv2["dv2_dfxdvz"] = self.abbre(sympy.diff(dv["dv_dfx"], self.vz))

        dv2["dv2_dfydX"] = self.abbre(sympy.diff(dv["dv_dfy"], self.X))
        dv2["dv2_dfydY"] = self.abbre(sympy.diff(dv["dv_dfy"], self.Y))
        dv2["dv2_dfydZ"] = self.abbre(sympy.diff(dv["dv_dfy"], self.Z))
        dv2["dv2_dfydfx"] = self.abbre(sympy.diff(dv["dv_dfy"], self.fx))
        dv2["dv2_dfydfy"] = self.abbre(sympy.diff(dv["dv_dfy"], self.fy))
        dv2["dv2_dfydcx"] = self.abbre(sympy.diff(dv["dv_dfy"], self.cx))
        dv2["dv2_dfydcy"] = self.abbre(sympy.diff(dv["dv_dfy"], self.cy))
        dv2["dv2_dfyds"] = self.abbre(sympy.diff(dv["dv_dfy"], self.s))
        dv2["dv2_dfydtx"] = self.abbre(sympy.diff(dv["dv_dfy"], self.tx))
        dv2["dv2_dfydty"] = self.abbre(sympy.diff(dv["dv_dfy"], self.ty))
        dv2["dv2_dfydtz"] = self.abbre(sympy.diff(dv["dv_dfy"], self.tz))
        dv2["dv2_dfydvx"] = self.abbre(sympy.diff(dv["dv_dfy"], self.vx))
        dv2["dv2_dfydvy"] = self.abbre(sympy.diff(dv["dv_dfy"], self.vy))
        dv2["dv2_dfydvz"] = self.abbre(sympy.diff(dv["dv_dfy"], self.vz))

        dv2["dv2_dcxdX"] = self.abbre(sympy.diff(dv["dv_dcx"], self.X))
        dv2["dv2_dcxdY"] = self.abbre(sympy.diff(dv["dv_dcx"], self.Y))
        dv2["dv2_dcxdZ"] = self.abbre(sympy.diff(dv["dv_dcx"], self.Z))
        dv2["dv2_dcxdfx"] = self.abbre(sympy.diff(dv["dv_dcx"], self.fx))
        dv2["dv2_dcxdfy"] = self.abbre(sympy.diff(dv["dv_dcx"], self.fy))
        dv2["dv2_dcxdcx"] = self.abbre(sympy.diff(dv["dv_dcx"], self.cx))
        dv2["dv2_dcxdcy"] = self.abbre(sympy.diff(dv["dv_dcx"], self.cy))
        dv2["dv2_dcxds"] = self.abbre(sympy.diff(dv["dv_dcx"], self.s))
        dv2["dv2_dcxdtx"] = self.abbre(sympy.diff(dv["dv_dcx"], self.tx))
        dv2["dv2_dcxdty"] = self.abbre(sympy.diff(dv["dv_dcx"], self.ty))
        dv2["dv2_dcxdtz"] = self.abbre(sympy.diff(dv["dv_dcx"], self.tz))
        dv2["dv2_dcxdvx"] = self.abbre(sympy.diff(dv["dv_dcx"], self.vx))
        dv2["dv2_dcxdvy"] = self.abbre(sympy.diff(dv["dv_dcx"], self.vy))
        dv2["dv2_dcxdvz"] = self.abbre(sympy.diff(dv["dv_dcx"], self.vz))

        dv2["dv2_dcydX"] = self.abbre(sympy.diff(dv["dv_dcy"], self.X))
        dv2["dv2_dcydY"] = self.abbre(sympy.diff(dv["dv_dcy"], self.Y))
        dv2["dv2_dcydZ"] = self.abbre(sympy.diff(dv["dv_dcy"], self.Z))
        dv2["dv2_dcydfx"] = self.abbre(sympy.diff(dv["dv_dcy"], self.fx))
        dv2["dv2_dcydfy"] = self.abbre(sympy.diff(dv["dv_dcy"], self.fy))
        dv2["dv2_dcydcx"] = self.abbre(sympy.diff(dv["dv_dcy"], self.cx))
        dv2["dv2_dcydcy"] = self.abbre(sympy.diff(dv["dv_dcy"], self.cy))
        dv2["dv2_dcyds"] = self.abbre(sympy.diff(dv["dv_dcy"], self.s))
        dv2["dv2_dcydtx"] = self.abbre(sympy.diff(dv["dv_dcy"], self.tx))
        dv2["dv2_dcydty"] = self.abbre(sympy.diff(dv["dv_dcy"], self.ty))
        dv2["dv2_dcydtz"] = self.abbre(sympy.diff(dv["dv_dcy"], self.tz))
        dv2["dv2_dcydvx"] = self.abbre(sympy.diff(dv["dv_dcy"], self.vx))
        dv2["dv2_dcydvy"] = self.abbre(sympy.diff(dv["dv_dcy"], self.vy))
        dv2["dv2_dcydvz"] = self.abbre(sympy.diff(dv["dv_dcy"], self.vz))

        dv2["dv2_dsdX"] = self.abbre(sympy.diff(dv["dv_ds"], self.X))
        dv2["dv2_dsdY"] = self.abbre(sympy.diff(dv["dv_ds"], self.Y))
        dv2["dv2_dsdZ"] = self.abbre(sympy.diff(dv["dv_ds"], self.Z))
        dv2["dv2_dsdfx"] = self.abbre(sympy.diff(dv["dv_ds"], self.fx))
        dv2["dv2_dsdfy"] = self.abbre(sympy.diff(dv["dv_ds"], self.fy))
        dv2["dv2_dsdcx"] = self.abbre(sympy.diff(dv["dv_ds"], self.cx))
        dv2["dv2_dsdcy"] = self.abbre(sympy.diff(dv["dv_ds"], self.cy))
        dv2["dv2_dsds"] = self.abbre(sympy.diff(dv["dv_ds"], self.s))
        dv2["dv2_dsdtx"] = self.abbre(sympy.diff(dv["dv_ds"], self.tx))
        dv2["dv2_dsdty"] = self.abbre(sympy.diff(dv["dv_ds"], self.ty))
        dv2["dv2_dsdtz"] = self.abbre(sympy.diff(dv["dv_ds"], self.tz))
        dv2["dv2_dsdvx"] = self.abbre(sympy.diff(dv["dv_ds"], self.vx))
        dv2["dv2_dsdvy"] = self.abbre(sympy.diff(dv["dv_ds"], self.vy))
        dv2["dv2_dsdvz"] = self.abbre(sympy.diff(dv["dv_ds"], self.vz))

        dv2["dv2_dtxdX"] = self.abbre(sympy.diff(dv["dv_dtx"], self.X))
        dv2["dv2_dtxdY"] = self.abbre(sympy.diff(dv["dv_dtx"], self.Y))
        dv2["dv2_dtxdZ"] = self.abbre(sympy.diff(dv["dv_dtx"], self.Z))
        dv2["dv2_dtxdfx"] = self.abbre(sympy.diff(dv["dv_dtx"], self.fx))
        dv2["dv2_dtxdfy"] = self.abbre(sympy.diff(dv["dv_dtx"], self.fy))
        dv2["dv2_dtxdcx"] = self.abbre(sympy.diff(dv["dv_dtx"], self.cx))
        dv2["dv2_dtxdcy"] = self.abbre(sympy.diff(dv["dv_dtx"], self.cy))
        dv2["dv2_dtxds"] = self.abbre(sympy.diff(dv["dv_dtx"], self.s))
        dv2["dv2_dtxdtx"] = self.abbre(sympy.diff(dv["dv_dtx"], self.tx))
        dv2["dv2_dtxdty"] = self.abbre(sympy.diff(dv["dv_dtx"], self.ty))
        dv2["dv2_dtxdtz"] = self.abbre(sympy.diff(dv["dv_dtx"], self.tz))
        dv2["dv2_dtxdvx"] = self.abbre(sympy.diff(dv["dv_dtx"], self.vx))
        dv2["dv2_dtxdvy"] = self.abbre(sympy.diff(dv["dv_dtx"], self.vy))
        dv2["dv2_dtxdvz"] = self.abbre(sympy.diff(dv["dv_dtx"], self.vz))

        dv2["dv2_dtydX"] = self.abbre(sympy.diff(dv["dv_dty"], self.X))
        dv2["dv2_dtydY"] = self.abbre(sympy.diff(dv["dv_dty"], self.Y))
        dv2["dv2_dtydZ"] = self.abbre(sympy.diff(dv["dv_dty"], self.Z))
        dv2["dv2_dtydfx"] = self.abbre(sympy.diff(dv["dv_dty"], self.fx))
        dv2["dv2_dtydfy"] = self.abbre(sympy.diff(dv["dv_dty"], self.fy))
        dv2["dv2_dtydcx"] = self.abbre(sympy.diff(dv["dv_dty"], self.cx))
        dv2["dv2_dtydcy"] = self.abbre(sympy.diff(dv["dv_dty"], self.cy))
        dv2["dv2_dtyds"] = self.abbre(sympy.diff(dv["dv_dty"], self.s))
        dv2["dv2_dtydtx"] = self.abbre(sympy.diff(dv["dv_dty"], self.tx))
        dv2["dv2_dtydty"] = self.abbre(sympy.diff(dv["dv_dty"], self.ty))
        dv2["dv2_dtydtz"] = self.abbre(sympy.diff(dv["dv_dty"], self.tz))
        dv2["dv2_dtydvx"] = self.abbre(sympy.diff(dv["dv_dty"], self.vx))
        dv2["dv2_dtydvy"] = self.abbre(sympy.diff(dv["dv_dty"], self.vy))
        dv2["dv2_dtydvz"] = self.abbre(sympy.diff(dv["dv_dty"], self.vz))

        dv2["dv2_dtzdX"] = self.abbre(sympy.diff(dv["dv_dtz"], self.X))
        dv2["dv2_dtzdY"] = self.abbre(sympy.diff(dv["dv_dtz"], self.Y))
        dv2["dv2_dtzdZ"] = self.abbre(sympy.diff(dv["dv_dtz"], self.Z))
        dv2["dv2_dtzdfx"] = self.abbre(sympy.diff(dv["dv_dtz"], self.fx))
        dv2["dv2_dtzdfy"] = self.abbre(sympy.diff(dv["dv_dtz"], self.fy))
        dv2["dv2_dtzdcx"] = self.abbre(sympy.diff(dv["dv_dtz"], self.cx))
        dv2["dv2_dtzdcy"] = self.abbre(sympy.diff(dv["dv_dtz"], self.cy))
        dv2["dv2_dtzds"] = self.abbre(sympy.diff(dv["dv_dtz"], self.s))
        dv2["dv2_dtzdtx"] = self.abbre(sympy.diff(dv["dv_dtz"], self.tx))
        dv2["dv2_dtzdty"] = self.abbre(sympy.diff(dv["dv_dtz"], self.ty))
        dv2["dv2_dtzdtz"] = self.abbre(sympy.diff(dv["dv_dtz"], self.tz))
        dv2["dv2_dtzdvx"] = self.abbre(sympy.diff(dv["dv_dtz"], self.vx))
        dv2["dv2_dtzdvy"] = self.abbre(sympy.diff(dv["dv_dtz"], self.vy))
        dv2["dv2_dtzdvz"] = self.abbre(sympy.diff(dv["dv_dtz"], self.vz))

        dv2["dv2_dvxdX"] = self.abbre(sympy.diff(dv["dv_dvx"], self.X))
        dv2["dv2_dvxdY"] = self.abbre(sympy.diff(dv["dv_dvx"], self.Y))
        dv2["dv2_dvxdZ"] = self.abbre(sympy.diff(dv["dv_dvx"], self.Z))
        dv2["dv2_dvxdfx"] = self.abbre(sympy.diff(dv["dv_dvx"], self.fx))
        dv2["dv2_dvxdfy"] = self.abbre(sympy.diff(dv["dv_dvx"], self.fy))
        dv2["dv2_dvxdcx"] = self.abbre(sympy.diff(dv["dv_dvx"], self.cx))
        dv2["dv2_dvxdcy"] = self.abbre(sympy.diff(dv["dv_dvx"], self.cy))
        dv2["dv2_dvxds"] = self.abbre(sympy.diff(dv["dv_dvx"], self.s))
        dv2["dv2_dvxdtx"] = self.abbre(sympy.diff(dv["dv_dvx"], self.tx))
        dv2["dv2_dvxdty"] = self.abbre(sympy.diff(dv["dv_dvx"], self.ty))
        dv2["dv2_dvxdtz"] = self.abbre(sympy.diff(dv["dv_dvx"], self.tz))
        dv2["dv2_dvxdvx"] = self.abbre(sympy.diff(dv["dv_dvx"], self.vx))
        dv2["dv2_dvxdvy"] = self.abbre(sympy.diff(dv["dv_dvx"], self.vy))
        dv2["dv2_dvxdvz"] = self.abbre(sympy.diff(dv["dv_dvx"], self.vz))

        dv2["dv2_dvydX"] = self.abbre(sympy.diff(dv["dv_dvy"], self.X))
        dv2["dv2_dvydY"] = self.abbre(sympy.diff(dv["dv_dvy"], self.Y))
        dv2["dv2_dvydZ"] = self.abbre(sympy.diff(dv["dv_dvy"], self.Z))
        dv2["dv2_dvydfx"] = self.abbre(sympy.diff(dv["dv_dvy"], self.fx))
        dv2["dv2_dvydfy"] = self.abbre(sympy.diff(dv["dv_dvy"], self.fy))
        dv2["dv2_dvydcx"] = self.abbre(sympy.diff(dv["dv_dvy"], self.cx))
        dv2["dv2_dvydcy"] = self.abbre(sympy.diff(dv["dv_dvy"], self.cy))
        dv2["dv2_dvyds"] = self.abbre(sympy.diff(dv["dv_dvy"], self.s))
        dv2["dv2_dvydtx"] = self.abbre(sympy.diff(dv["dv_dvy"], self.tx))
        dv2["dv2_dvydty"] = self.abbre(sympy.diff(dv["dv_dvy"], self.ty))
        dv2["dv2_dvydtz"] = self.abbre(sympy.diff(dv["dv_dvy"], self.tz))
        dv2["dv2_dvydvx"] = self.abbre(sympy.diff(dv["dv_dvy"], self.vx))
        dv2["dv2_dvydvy"] = self.abbre(sympy.diff(dv["dv_dvy"], self.vy))
        dv2["dv2_dvydvz"] = self.abbre(sympy.diff(dv["dv_dvy"], self.vz))

        dv2["dv2_dvzdX"] = self.abbre(sympy.diff(dv["dv_dvz"], self.X))
        dv2["dv2_dvzdY"] = self.abbre(sympy.diff(dv["dv_dvz"], self.Y))
        dv2["dv2_dvzdZ"] = self.abbre(sympy.diff(dv["dv_dvz"], self.Z))
        dv2["dv2_dvzdfx"] = self.abbre(sympy.diff(dv["dv_dvz"], self.fx))
        dv2["dv2_dvzdfy"] = self.abbre(sympy.diff(dv["dv_dvz"], self.fy))
        dv2["dv2_dvzdcx"] = self.abbre(sympy.diff(dv["dv_dvz"], self.cx))
        dv2["dv2_dvzdcy"] = self.abbre(sympy.diff(dv["dv_dvz"], self.cy))
        dv2["dv2_dvzds"] = self.abbre(sympy.diff(dv["dv_dvz"], self.s))
        dv2["dv2_dvzdtx"] = self.abbre(sympy.diff(dv["dv_dvz"], self.tx))
        dv2["dv2_dvzdty"] = self.abbre(sympy.diff(dv["dv_dvz"], self.ty))
        dv2["dv2_dvzdtz"] = self.abbre(sympy.diff(dv["dv_dvz"], self.tz))
        dv2["dv2_dvzdvx"] = self.abbre(sympy.diff(dv["dv_dvz"], self.vx))
        dv2["dv2_dvzdvy"] = self.abbre(sympy.diff(dv["dv_dvz"], self.vy))
        dv2["dv2_dvzdvz"] = self.abbre(sympy.diff(dv["dv_dvz"], self.vz))


        return du2, dv2

    def abbre_normal_angle(self, expr):
        return expr.subs(sympy.sqrt(self.vx**2 + self.vy**2 + self.vz**2), 'theta')

    def abbre_small_angle(self, expr):
        return expr.subs(sympy.sqrt(self.vx**2 + self.vy**2 + self.vz**2), 'theta')


