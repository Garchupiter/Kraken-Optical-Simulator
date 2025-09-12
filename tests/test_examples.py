import unittest
import sys
import subprocess


class TestExamples(unittest.TestCase):

    def test_Axicon(self):
        """Test Examp_Axicon.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_Axicon.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in Axicon.py: {result.stderr}")
        print(f"Examp_Axicon.py output:", result.stdout.strip())

    def test_Axicon_And_Cylinder(self):
        """Test Examp_Axicon_And_Cylinder.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_Axicon_And_Cylinder.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in Axicon_And_Cylinder.py: {result.stderr}")
        print(f"Examp_Axicon_And_Cylinder.py output:", result.stdout.strip())

    def test_CzernyTurner(self):
        """Test Examp_CzernyTurner.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_CzernyTurner.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in CzernyTurner.py: {result.stderr}")
        print(f"Examp_CzernyTurner.py output:", result.stdout.strip())

    def test_Perfect_lens(self):
        """Test Examp_Perfect_lens.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_Perfect_lens.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in Perfect_lens.py: {result.stderr}")
        print(f"Examp_Perfect_lens.py output:", result.stdout.strip())

    def test_Perfect_lens_Telescope(self):
        """Test Examp_Perfect_lens_Telescope.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_Perfect_lens_Telescope.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in Perfect_lens_Telescope.py: {result.stderr}")
        print(f"Examp_Perfect_lens_Telescope.py output:", result.stdout.strip())

    def test_Diffraction_Grating_Reflection(self):
        """Test Examp_Diffraction_Grating_Reflection.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_Diffraction_Grating_Reflection.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in Diffraction_Grating_Reflection.py: {result.stderr}")
        print(f"Examp_Diffraction_Grating_Reflection.py output:", result.stdout.strip())

    def test_Diffraction_Grating_Reflection_Single(self):
        """Test Examp_Diffraction_Grating_Reflection_Single.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_Diffraction_Grating_Reflection_Single.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in Diffraction_Grating_Reflection_Single.py: {result.stderr}")
        print(f"Examp_Diffraction_Grating_Reflection_Single.py output:", result.stdout.strip())

    def test_Diffraction_Grating_Transmission(self):
        """Test Examp_Diffraction_Grating_Transmission.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_Diffraction_Grating_Transmission.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in Diffraction_Grating_Transmission.py: {result.stderr}")
        print(f"Examp_Diffraction_Grating_Transmission.py output:", result.stdout.strip())

    def test_Dispersion_By_AbbeNumber(self):
        """Test Examp_Dispersion_By_AbbeNumber.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_Dispersion_By_AbbeNumber.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in Dispersion_By_AbbeNumber.py: {result.stderr}")
        print(f"Examp_Dispersion_By_AbbeNumber.py output:", result.stdout.strip())

    def test_Doublet_Lens_ParaxMatrix(self):
        """Test Examp_Doublet_Lens-ParaxMatrix.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_Doublet_Lens-ParaxMatrix.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in Doublet_Lens-ParaxMatrix.py: {result.stderr}")
        print(f"Examp_Doublet_Lens-ParaxMatrix.py output:", result.stdout.strip())

    def test_Doublet_Lens(self):
        """Test Examp_Doublet_Lens.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_Doublet_Lens.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in Doublet_Lens.py: {result.stderr}")
        print(f"Examp_Doublet_Lens.py output:", result.stdout.strip())

    def test_Doublet_Lens_3Dcolor(self):
        """Test Examp_Doublet_Lens_3Dcolor.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_Doublet_Lens_3Dcolor.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in Doublet_Lens_3Dcolor.py: {result.stderr}")
        print(f"Examp_Doublet_Lens_3Dcolor.py output:", result.stdout.strip())

    def test_Doublet_Lens_CommandsSystem(self):
        """Test Examp_Doublet_Lens_CommandsSystem.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_Doublet_Lens_CommandsSystem.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in Doublet_Lens_CommandsSystem.py: {result.stderr}")
        print(f"Examp_Doublet_Lens_CommandsSystem.py output:", result.stdout.strip())

    def test_Doublet_Lens_Cylinder(self):
        """Test Examp_Doublet_Lens_Cylinder.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_Doublet_Lens_Cylinder.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in Doublet_Lens_Cylinder.py: {result.stderr}")
        print(f"Examp_Doublet_Lens_Cylinder.py output:", result.stdout.strip())

    def test_Doublet_Lens_NonSec_AR_Coating(self):
        """Test Examp_Doublet_Lens_NonSec-AR_Coating.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_Doublet_Lens_NonSec-AR_Coating.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in Doublet_Lens_NonSec-AR_Coating.py: {result.stderr}")
        print(f"Examp_Doublet_Lens_NonSec-AR_Coating.py output:", result.stdout.strip())

    def test_Doublet_Lens_NonSec(self):
        """Test Examp_Doublet_Lens_NonSec.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_Doublet_Lens_NonSec.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in Doublet_Lens_NonSec.py: {result.stderr}")
        print(f"Examp_Doublet_Lens_NonSec.py output:", result.stdout.strip())

    def test_Doublet_Lens_Pupil(self):
        """Test Examp_Doublet_Lens_Pupil.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_Doublet_Lens_Pupil.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in Doublet_Lens_Pupil.py: {result.stderr}")
        print(f"Examp_Doublet_Lens_Pupil.py output:", result.stdout.strip())

    def test_Doublet_Lens_Pupil_Seidel(self):
        """Test Examp_Doublet_Lens_Pupil_Seidel.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_Doublet_Lens_Pupil_Seidel.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in Doublet_Lens_Pupil_Seidel.py: {result.stderr}")
        print(f"Examp_Doublet_Lens_Pupil_Seidel.py output:", result.stdout.strip())

    def test_Doublet_Lens_Tilt_Nulls(self):
        """Test Examp_Doublet_Lens_Tilt-Nulls.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_Doublet_Lens_Tilt-Nulls.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in Doublet_Lens_Tilt-Nulls.py: {result.stderr}")
        print(f"Examp_Doublet_Lens_Tilt-Nulls.py output:", result.stdout.strip())

    def test_Doublet_Lens_Tilt(self):
        """Test Examp_Doublet_Lens_Tilt.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_Doublet_Lens_Tilt.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in Doublet_Lens_Tilt.py: {result.stderr}")
        print(f"Examp_Doublet_Lens_Tilt.py output:", result.stdout.strip())

    def test_Doublet_Lens_Tilt_non_sec_AR_Coating(self):
        """Test Examp_Doublet_Lens_Tilt_non_sec-AR-Coating.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_Doublet_Lens_Tilt_non_sec-AR-Coating.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in Doublet_Lens_Tilt_non_sec-AR-Coating.py: {result.stderr}")
        print(f"Examp_Doublet_Lens_Tilt_non_sec-AR-Coating.py output:", result.stdout.strip())

    def test_Doublet_Lens_Tilt_non_sec(self):
        """Test Examp_Doublet_Lens_Tilt_non_sec.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_Doublet_Lens_Tilt_non_sec.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in Doublet_Lens_Tilt_non_sec.py: {result.stderr}")
        print(f"Examp_Doublet_Lens_Tilt_non_sec.py output:", result.stdout.strip())

    def test_Doublet_Lens_Zernike(self):
        """Test Examp_Doublet_Lens_Zernike.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_Doublet_Lens_Zernike.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in Doublet_Lens_Zernike.py: {result.stderr}")
        print(f"Examp_Doublet_Lens_Zernike.py output:", result.stdout.strip())

    def test_Doublet_Optimization(self):
        """Test Examp_Doublet_Optimization.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_Doublet_Optimization.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in Doublet_Optimization.py: {result.stderr}")
        print(f"Examp_Doublet_Optimization.py output:", result.stdout.strip())

    def test_ExtraShape_Micro_Lens_Array(self):
        """Test Examp_ExtraShape_Micro_Lens_Array.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_ExtraShape_Micro_Lens_Array.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in ExtraShape_Micro_Lens_Array.py: {result.stderr}")
        print(f"Examp_ExtraShape_Micro_Lens_Array.py output:", result.stdout.strip())

    def test_ExtraShape_Radial_Sine(self):
        """Test Examp_ExtraShape_Radial_Sine.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_ExtraShape_Radial_Sine.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in ExtraShape_Radial_Sine.py: {result.stderr}")
        print(f"Examp_ExtraShape_Radial_Sine.py output:", result.stdout.strip())

    def test_ExtraShape_UserFacets(self):
        """Test Examp_ExtraShape_UserFacets.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_ExtraShape_UserFacets.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in ExtraShape_UserFacets.py: {result.stderr}")
        print(f"Examp_ExtraShape_UserFacets.py output:", result.stdout.strip())

    def test_ExtraShape_XY_Cosines(self):
        """Test Examp_ExtraShape_XY_Cosines.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_ExtraShape_XY_Cosines.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in ExtraShape_XY_Cosines.py: {result.stderr}")
        print(f"Examp_ExtraShape_XY_Cosines.py output:", result.stdout.strip())

    def test_ExtraShape_XY_Cosines_UDA(self):
        """Test Examp_ExtraShape_XY_Cosines_UDA.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_ExtraShape_XY_Cosines_UDA.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in ExtraShape_XY_Cosines_UDA.py: {result.stderr}")
        print(f"Examp_ExtraShape_XY_Cosines_UDA.py output:", result.stdout.strip())

    def test_Flat_Mirror_45Deg(self):
        """Test Examp_Flat_Mirror_45Deg.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_Flat_Mirror_45Deg.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in Flat_Mirror_45Deg.py: {result.stderr}")
        print(f"Examp_Flat_Mirror_45Deg.py output:", result.stdout.strip())

    def test_Flat_NonSec_AR_caoating(self):
        """Test Examp_Flat_NonSec_AR-caoating.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_Flat_NonSec_AR-caoating.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in Flat_NonSec_AR-caoating.py: {result.stderr}")
        print(f"Examp_Flat_NonSec_AR-caoating.py output:", result.stdout.strip())

    def test_Fresnel(self):
        """Test Examp_Fresnel.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_Fresnel.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in Fresnel.py: {result.stderr}")
        print(f"Examp_Fresnel.py output:", result.stdout.strip())

    def test_ParaboleMirror_Shift(self):
        """Test Examp_ParaboleMirror_Shift.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_ParaboleMirror_Shift.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in ParaboleMirror_Shift.py: {result.stderr}")
        print(f"Examp_ParaboleMirror_Shift.py output:", result.stdout.strip())

    def test_ParaboleMirror_Shift_UDA(self):
        """Test Examp_ParaboleMirror_Shift_UDA.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_ParaboleMirror_Shift_UDA.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in ParaboleMirror_Shift_UDA.py: {result.stderr}")
        print(f"Examp_ParaboleMirror_Shift_UDA.py output:", result.stdout.strip())


    def test_Pickle_Doublet_Lens_3Dcolor(self):
        """Test Examp_Pickle_Doublet_Lens_3Dcolor.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_Pickle_Doublet_Lens_3Dcolor.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in Pickle_Doublet_Lens_3Dcolor.py: {result.stderr}")
        print(f"Examp_Pickle_Doublet_Lens_3Dcolor.py output:", result.stdout.strip())

    def test_Prism_STL_AR_coating(self):
        """Test Examp_Prism_STL-AR_coating.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_Prism_STL-AR_coating.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in Prism_STL-AR_coating.py: {result.stderr}")
        print(f"Examp_Prism_STL-AR_coating.py output:", result.stdout.strip())

    def test_Prism_STL(self):
        """Test Examp_Prism_STL.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_Prism_STL.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in Prism_STL.py: {result.stderr}")
        print(f"Examp_Prism_STL.py output:", result.stdout.strip())

    def test_Ray(self):
        """Test Examp_Ray.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_Ray.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in Ray.py: {result.stderr}")
        print(f"Examp_Ray.py output:", result.stdout.strip())

    def test_Refraction_Prism(self):
        """Test Examp_Refraction_Prism.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_Refraction_Prism.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in Refraction_Prism.py: {result.stderr}")
        print(f"Examp_Refraction_Prism.py output:", result.stdout.strip())

    def test_Refraction_Prism_OneRay(self):
        """Test Examp_Refraction_Prism_OneRay.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_Refraction_Prism_OneRay.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in Refraction_Prism_OneRay.py: {result.stderr}")
        print(f"Examp_Refraction_Prism_OneRay.py output:", result.stdout.strip())

    def test_Refraction_Prism_solid(self):
        """Test Examp_Refraction_Prism_solid.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_Refraction_Prism_solid.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in Refraction_Prism_solid.py: {result.stderr}")
        print(f"Examp_Refraction_Prism_solid.py output:", result.stdout.strip())

    def test_Refraction_Prism_solid_Generation(self):
        """Test Examp_Refraction_Prism_solid_Generation.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_Refraction_Prism_solid_Generation.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in Refraction_Prism_solid_Generation.py: {result.stderr}")
        print(f"Examp_Refraction_Prism_solid_Generation.py output:", result.stdout.strip())

    def test_RonchiTest(self):
        """Test Examp_RonchiTest.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_RonchiTest.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in RonchiTest.py: {result.stderr}")
        print(f"Examp_RonchiTest.py output:", result.stdout.strip())

    def test_Solid_Object_STL(self):
        """Test Examp_Solid_Object_STL.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_Solid_Object_STL.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in Solid_Object_STL.py: {result.stderr}")
        print(f"Examp_Solid_Object_STL.py output:", result.stdout.strip())

    def test_Solid_Object_STL_ARRAY(self):
        """Test Examp_Solid_Object_STL_ARRAY.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_Solid_Object_STL_ARRAY.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in Solid_Object_STL_ARRAY.py: {result.stderr}")
        print(f"Examp_Solid_Object_STL_ARRAY.py output:", result.stdout.strip())

    def test_Source_Distribution_Function(self):
        """Test Examp_Source_Distribution_Function.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_Source_Distribution_Function.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in Source_Distribution_Function.py: {result.stderr}")
        print(f"Examp_Source_Distribution_Function.py output:", result.stdout.strip())

    def test_Sphere_2(self):
        """Test Examp_Sphere 2.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_Sphere 2.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in Sphere 2.py: {result.stderr}")
        print(f"Examp_Sphere 2.py output:", result.stdout.strip())

    def test_Sphere(self):
        """Test Examp_Sphere.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_Sphere.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in Sphere.py: {result.stderr}")
        print(f"Examp_Sphere.py output:", result.stdout.strip())

    def test_Spruce_tone_Github_User_Manually_enter_the_refractive_index_dispersion_and_alpha(self):
        """Test Examp_Spruce-tone_Github_User (Manually enter the refractive index, dispersion and alpha).py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_Spruce-tone_Github_User (Manually enter the refractive index, dispersion and alpha).py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in Spruce-tone_Github_User (Manually enter the refractive index, dispersion and alpha).py: {result.stderr}")
        print(f"Examp_Spruce-tone_Github_User (Manually enter the refractive index, dispersion and alpha).py output:", result.stdout.strip())

    def test_Spruce_tone_Github_User_Loading_Zemax_and_Catalogs(self):
        """Test Examp_Spruce-tone_Github_User(Loading_Zemax_and_Catalogs).py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_Spruce-tone_Github_User(Loading_Zemax_and_Catalogs).py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in Spruce-tone_Github_User(Loading_Zemax_and_Catalogs).py: {result.stderr}")
        print(f"Examp_Spruce-tone_Github_User(Loading_Zemax_and_Catalogs).py output:", result.stdout.strip())

    def test_Tel_2M_STL_ImageSlicer(self):
        """Test Examp_Tel_2M-STL_ImageSlicer.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_Tel_2M-STL_ImageSlicer.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in Tel_2M-STL_ImageSlicer.py: {result.stderr}")
        print(f"Examp_Tel_2M-STL_ImageSlicer.py output:", result.stdout.strip())

    def test_Tel_2M(self):
        """Test Examp_Tel_2M.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_Tel_2M.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in Tel_2M.py: {result.stderr}")
        print(f"Examp_Tel_2M.py output:", result.stdout.strip())

    def test_Tel_2M_Atmospheric_Refraction_Corrector_Adaptable(self):
        """Test Examp_Tel_2M_Atmospheric_Refraction_Corrector_Adaptable.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_Tel_2M_Atmospheric_Refraction_Corrector_Adaptable.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in Tel_2M_Atmospheric_Refraction_Corrector_Adaptable.py: {result.stderr}")
        print(f"Examp_Tel_2M_Atmospheric_Refraction_Corrector_Adaptable.py output:", result.stdout.strip())

    def test_Tel_2M_Atmospheric_Refraction_Corrector_Static(self):
        """Test Examp_Tel_2M_Atmospheric_Refraction_Corrector_Static.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_Tel_2M_Atmospheric_Refraction_Corrector_Static.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in Tel_2M_Atmospheric_Refraction_Corrector_Static.py: {result.stderr}")
        print(f"Examp_Tel_2M_Atmospheric_Refraction_Corrector_Static.py output:", result.stdout.strip())

    def test_Tel_2M_Cuña(self):
        """Test Examp_Tel_2M_Cuña.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_Tel_2M_Cuña.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in Tel_2M_Cuña.py: {result.stderr}")
        print(f"Examp_Tel_2M_Cuña.py output:", result.stdout.strip())

    def test_Tel_2M_Echelle(self):
        """Test Examp_Tel_2M_Echelle.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_Tel_2M_Echelle.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in Tel_2M_Echelle.py: {result.stderr}")
        print(f"Examp_Tel_2M_Echelle.py output:", result.stdout.strip())

    def test_Tel_2M_Optimization_Variables(self):
        """Test Examp_Tel_2M_Optimization_Variables.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_Tel_2M_Optimization_Variables.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in Tel_2M_Optimization_Variables.py: {result.stderr}")
        print(f"Examp_Tel_2M_Optimization_Variables.py output:", result.stdout.strip())

    def test_Tel_2M_Pupila(self):
        """Test Examp_Tel_2M_Pupila.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_Tel_2M_Pupila.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in Tel_2M_Pupila.py: {result.stderr}")
        print(f"Examp_Tel_2M_Pupila.py output:", result.stdout.strip())

    def test_Tel_2M_Spyder_Spot_Diagram(self):
        """Test Examp_Tel_2M_Spyder_Spot_Diagram.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_Tel_2M_Spyder_Spot_Diagram.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in Tel_2M_Spyder_Spot_Diagram.py: {result.stderr}")
        print(f"Examp_Tel_2M_Spyder_Spot_Diagram.py output:", result.stdout.strip())

    def test_Tel_2M_Spyder_Spot_RMS(self):
        """Test Examp_Tel_2M_Spyder_Spot_RMS.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_Tel_2M_Spyder_Spot_RMS.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in Tel_2M_Spyder_Spot_RMS.py: {result.stderr}")
        print(f"Examp_Tel_2M_Spyder_Spot_RMS.py output:", result.stdout.strip())

    def test_Tel_2M_Wavefront_Fitting(self):
        """Test Examp_Tel_2M_Wavefront_Fitting.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_Tel_2M_Wavefront_Fitting.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in Tel_2M_Wavefront_Fitting.py: {result.stderr}")
        print(f"Examp_Tel_2M_Wavefront_Fitting.py output:", result.stdout.strip())

    def test_Tel_2M_Wavefront_Fitting_optimization(self):
        """Test Examp_Tel_2M_Wavefront_Fitting_optimization.py"""
        result = subprocess.run([sys.executable, "KrakenOS/Examples/Examp_Tel_2M_Wavefront_Fitting_optimization.py"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, msg=f"Error in Tel_2M_Wavefront_Fitting_optimization.py: {result.stderr}")
        print(f"Examp_Tel_2M_Wavefront_Fitting_optimization.py output:", result.stdout.strip())
