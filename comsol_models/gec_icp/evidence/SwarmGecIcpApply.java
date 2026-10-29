import com.comsol.model.*;
import com.comsol.model.util.*;

public class SwarmGecIcpApply {
  public static void main(String[] args) throws Exception {
    Model model = ModelUtil.loadCopy("swarmGecIcpApply", "C:/Users/user/Desktop/swarm_5/comsol_modes/argon_gec_icp.mph");
    // Restricted local-mean-energy closure shared by all
    // Swarm layers: external mu*N and f0(E,<E>); COMSOL
    // retains diffusion, energy transport, RF and chemistry.
    model.component("comp1").physics("plas").prop("ElectronProperties").set("ReducedProps", true);
    model.component("comp1").physics("plas").prop("ElectronProperties").set("TensorElectronProps", false);
    model.component("comp1").physics("plas").prop("ElectronProperties").set("IncludeThermalDiffusion", false);
    model.component("comp1").physics("plas").prop("ElectronProperties").set("MeanElectronEnergyModel", "LocalEnergyApproximationE");
    model.component("comp1").physics("plas").feature("pes1").set("SpecifyElectronDensityAndEnergy", "SpecifyMueOnly");
    try { model.component("comp1").variable("swIcpClosureVars"); } catch (Exception missing) {
      model.component("comp1").variable().create("swIcpClosureVars");
    }
    model.component("comp1").variable("swIcpClosureVars").set("sw_icp_logeps", "(En-Ne)");
    replaceInterpolation(model, "sw_icp_log_muN");
    model.func("sw_icp_log_muN").set("source", "table");
    model.func("sw_icp_log_muN").set("nargs", 1);
    model.func("sw_icp_log_muN").set("argunit", "1");
    model.func("sw_icp_log_muN").set("fununit", "1");
    model.func("sw_icp_log_muN").set("interp", "piecewisecubic");
    model.func("sw_icp_log_muN").set("extrap", "const");
    model.func("sw_icp_log_muN").set("table", new String[][]{new String[]{"8.99344625992673397e-01", "5.63673893112150353e+01"}, new String[]{"1.24176720657430084e+00", "5.58460227541406553e+01"}, new String[]{"1.43351048362289846e+00", "5.55349520672534922e+01"}, new String[]{"1.61108962922305721e+00", "5.53023549793510156e+01"}, new String[]{"1.67869359068208168e+00", "5.53233412942852283e+01"}, new String[]{"1.72854186000140730e+00", "5.52742178382623379e+01"}, new String[]{"1.76882367672341334e+00", "5.52170689133086938e+01"}, new String[]{"1.83203660433106319e+00", "5.51382539941497498e+01"}, new String[]{"1.93606221753196039e+00", "5.50462721644598858e+01"}, new String[]{"2.16675546069127600e+00", "5.49200739213975879e+01"}});
    model.component("comp1").physics("plas").feature("pes1").set("muN", new String[]{"exp(sw_icp_log_muN(sw_icp_logeps))*1[1/(V*m*s)]", "0[1/(V*m*s)]", "0[1/(V*m*s)]", "0[1/(V*m*s)]", "exp(sw_icp_log_muN(sw_icp_logeps))*1[1/(V*m*s)]", "0[1/(V*m*s)]", "0[1/(V*m*s)]", "0[1/(V*m*s)]", "exp(sw_icp_log_muN(sw_icp_logeps))*1[1/(V*m*s)]"});
    try { model.component("comp1").func().remove("sw_icp_eedf_two_term"); } catch (Exception missing) {}
    model.component("comp1").func().create("sw_icp_eedf_two_term", "Interpolation");
    model.component("comp1").func("sw_icp_eedf_two_term").set("source", "file");
    model.component("comp1").func("sw_icp_eedf_two_term").set("filename", "C:/Users/user/Desktop/swarm_5/comsol_models/gec_icp/bundle/eedf_f0_comsol_2d.csv");
    model.component("comp1").func("sw_icp_eedf_two_term").set("nargs", 2);
    model.component("comp1").func("sw_icp_eedf_two_term").set("struct", "spreadsheet");
    model.component("comp1").func("sw_icp_eedf_two_term").set("scaledata", "auto");
    model.component("comp1").func("sw_icp_eedf_two_term").set("funcnametable", new String[][]{new String[]{"sw_icp_eedf_two_term", "1"}});
    model.component("comp1").func("sw_icp_eedf_two_term").set("argunit", "eV,eV");
    model.component("comp1").func("sw_icp_eedf_two_term").set("fununit", "1");
    model.component("comp1").func("sw_icp_eedf_two_term").set("interp", "linear");
    model.component("comp1").func("sw_icp_eedf_two_term").set("extrap", "const");
    model.component("comp1").func("sw_icp_eedf_two_term").importData();
    model.component("comp1").physics("plas").prop("EEDFSettings").set("eedf", "sw_icp_eedf_two_term");
    model.component("comp1").physics("plas").feature("eir1").set("SpecifyReactionUsing", "UseCrossSectionData");
    model.component("comp1").physics("plas").feature("eir1").set("eedf", "FromPhysicsInterfaceProperty");
    model.component("comp1").physics("plas").feature("eir1").set("UseTownsendTwoTermsBoltzmann", false);
    model.component("comp1").physics("plas").feature("eir2").set("SpecifyReactionUsing", "UseCrossSectionData");
    model.component("comp1").physics("plas").feature("eir2").set("eedf", "FromPhysicsInterfaceProperty");
    model.component("comp1").physics("plas").feature("eir2").set("UseTownsendTwoTermsBoltzmann", false);
    model.component("comp1").physics("plas").feature("eir3").set("SpecifyReactionUsing", "UseCrossSectionData");
    model.component("comp1").physics("plas").feature("eir3").set("eedf", "FromPhysicsInterfaceProperty");
    model.component("comp1").physics("plas").feature("eir3").set("UseTownsendTwoTermsBoltzmann", false);
    model.component("comp1").physics("plas").feature("eir4").set("SpecifyReactionUsing", "UseCrossSectionData");
    model.component("comp1").physics("plas").feature("eir4").set("eedf", "FromPhysicsInterfaceProperty");
    model.component("comp1").physics("plas").feature("eir4").set("UseTownsendTwoTermsBoltzmann", false);
    model.component("comp1").physics("plas").feature("eir5").set("SpecifyReactionUsing", "UseCrossSectionData");
    model.component("comp1").physics("plas").feature("eir5").set("eedf", "FromPhysicsInterfaceProperty");
    model.component("comp1").physics("plas").feature("eir5").set("UseTownsendTwoTermsBoltzmann", false);
    model.save("C:/Users/user/Desktop/swarm_5/outputs/configured_models/gec_icp/solved.mph");
    ModelUtil.remove("swarmGecIcpApply");
  }

  private static void replaceInterpolation(Model model, String tag) {
    try { model.func().remove(tag); } catch (Exception ignored) {}
    model.func().create(tag, "Interpolation");
  }
}
