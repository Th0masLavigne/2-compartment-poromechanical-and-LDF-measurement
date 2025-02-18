from .constitutive_laws import readme_CL
from .variational_forms import readme_VF
from .variational_forms_Neo_Hooke import readme_VFNH
from .variational_forms_TW import readme_VFTW
from .functions import readme_F

__all__ = ["constitutive_laws", "functions", "variational_forms",  "variational_forms_TW", "variational_forms_Neo_Hooke"]

if __name__=='__main__':
	print("The package 'Porous_media_package' was successfully loaded")