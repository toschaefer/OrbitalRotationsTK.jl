module OrbitalRotationsTK

using LinearAlgebra
import Random
using DFTK
import PsiTK
import Lucon


export Monomial
export CustomFunction
public derivative
public taylor_degree
include("common/scalar_functions.jl")

export OrbitalSubspace
export RealSpace
export FourierSpace
public prepare_gradient
public gradient
public maximize
public max_taylor_degree
include("common/representations.jl")

public JointDiagonalizationFunctional
public one_body_operators
include("functionals/joint_diagonalization/joint_diagonalization.jl")

export NPL
include("functionals/joint_diagonalization/npl.jl")

export RotationResult
export rotate
include("rotate.jl")


end
