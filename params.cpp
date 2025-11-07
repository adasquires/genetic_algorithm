#include "genetic_algorithm.h"

std::vector<std::string> param_names =
    {"alpha", "beta", "FoodPos_x", "foodPos_y", "gamma", "kappa", "lambda",
    "AWA_AIY", "AIY_AIY", "AIY_RIA", "RIA_RMDD", "RIA_RMDV", "SMDD_RIA", "SMDV_RIA", "RIA_RIA",
    "RIA_SMDD", "RIA_SMDV", "RMDD_RIA", "RMDV_RIA", "SMDD_SMDD", "SMDV_SMDV", "RMDD_RMDD", "RMDV_RMDV",
    "SMDD_SMDV", "SMDV_SMDD", "SMDD_RMDV", "SMDV_RMDD", "RMDD_RMDV", "RMDV_RMDD", "SMDD_RMDV_ele", "SMDV_RMDV_ele", "RMDV_RMDV_ele",
    "AIY_tau", "AIY_theta", "AWA_tau", "AWA_theta", "RIA_tau", "RIA_theta",
    "RMDD_tau", "RMDD_theta", "RMDV_tau", "RMDV_theta", "SMDD_tau", "SMDD_theta", "SMDV_tau", "SMDV_theta",
    "NMJ_DB", "NMJ_DD", "NMJ_RMDD", "NMJ_RMDV", "NMJ_SMDD", "NMJ_SMDV",
    "NMJ_VBA", "NMJ_VBP", "NMJ_VDA", "NMJ_VDP",
    "sr_headgain", "sr_vcgain",
    "DB_DB", "VBA_VBA", "VBP_VBP", "DD_DD", "VDA_VDA", "VDP_VDP", "DB_DD",
    "VBA_VDA", "VBP_VDP", "DB_VDA", "DB_VDP", "VBA_DD", "VBP_DD", "DD_VDA",
    "DD_VDA_ele", "DD_VDP_ele", "VDA_VDP_ele", "VBA_VBP_ele",
    "fwd_DB_DB", "fwd_VBP_VBA", "fwd_DD_DD", "fwd_VDP_VDA", "fwd_VBP_DB",
    "DB_tau", "DB_theta", "DD_tau", "DD_theta", "VBA_tau", "VBA_theta",
    "VBP_tau", "VBP_theta", "VDA_tau", "VDA_theta", "VDP_tau", "VDP_theta"};

std::vector<std::string> change_params =
    {"alpha", "beta", "foodPos_x", "foodPos_y", "gamma",
     "AWA_AIY", "AIY_AIY", "AIY_RIA", "RIA_RMDD", "RIA_RMDV", "SMDD_RIA", "SMDV_RIA", "RIA_RIA",
     "RIA_SMDD", "RIA_SMDV", "RMDD_RIA", "RMDV_RIA", "SMDD_SMDD", "SMDV_SMDV", "RMDD_RMDD", "RMDV_RMDV",
     "SMDD_SMDV", "SMDV_SMDD", "SMDD_RMDV", "SMDV_RMDD", "RMDD_RMDV", "RMDV_RMDD", "SMDD_RMDD_ele", "SMDV_RMDV_ele", "RMDV_RMDD_ele",
     "AIY_tau", "AIY_theta", "AWA_tau", "AWA_theta", "RIA_tau", "RIA_theta",
     "RMDD_tau", "RMDD_theta", "RMDV_tau", "RMDV_theta", "SMDD_tau", "SMDD_theta", "SMDV_tau", "SMDV_theta"};

// Parse parameters from input/params.json.
std::map<std::string, double> parse_params(json p, std::string input) {

    std::map<std::string, double> params;

    std::ifstream f(input);

    if (!f) {
        std::cerr << "params.json not found" << std::endl;
    }

    json data;
    f >> data;

    // Chemoreceptor parameters.
    params["alpha"] = data["ChemoReceptors"]["alpha"];
    params["beta"] = data["ChemoReceptors"]["beta"];
    params["foodPos_x"] = data["ChemoReceptors"]["foodPos"]["x"];
    params["foodPos_y"] = data["ChemoReceptors"]["foodPos"]["y"];
    params["gamma"] = data["ChemoReceptors"]["gamma"];
    params["kappa"] = data["ChemoReceptors"]["kappa"];
    params["lambda"] = data["ChemoReceptors"]["lambda"];

    // Head connection parameters.
    params["AWA_AIY"] = data["Head"]["connections"][0]["weight"];
    params["AIY_AIY"] = data["Head"]["connections"][1]["weight"];
    params["AIY_RIA"] = data["Head"]["connections"][2]["weight"];
    params["RIA_RMDD"] = data["Head"]["connections"][3]["weight"];
    params["RIA_RMDV"] = data["Head"]["connections"][4]["weight"];
    params["SMDD_RIA"] = data["Head"]["connections"][5]["weight"];
    params["SMDV_RIA"] = data["Head"]["connections"][6]["weight"];
    params["RIA_RIA"] = data["Head"]["connections"][7]["weight"];
    params["RIA_SMDD"] = data["Head"]["connections"][8]["weight"];
    params["RIA_SMDV"] = data["Head"]["connections"][9]["weight"];
    params["RMDD_RIA"] = data["Head"]["connections"][10]["weight"];
    params["RMDV_RIA"] = data["Head"]["connections"][11]["weight"];
    params["SMDD_SMDD"] = data["Head"]["connections"][12]["weight"];
    params["SMDV_SMDV"] = data["Head"]["connections"][13]["weight"];
    params["RMDD_RMDD"] = data["Head"]["connections"][14]["weight"];
    params["RMDV_RMDV"] = data["Head"]["connections"][15]["weight"];
    params["SMDD_SMDV"] = data["Head"]["connections"][16]["weight"];
    params["SMDV_SMDD"] = data["Head"]["connections"][17]["weight"];
    params["SMDD_RMDV"] = data["Head"]["connections"][18]["weight"];
    params["SMDV_RMDD"] = data["Head"]["connections"][19]["weight"];
    params["RMDD_RMDV"] = data["Head"]["connections"][20]["weight"];
    params["RMDV_RMDD"] = data["Head"]["connections"][21]["weight"];
    params["SMDD_RMDD_ele"] = data["Head"]["connections"][22]["weight"];
    params["SMDV_RMDV_ele"] = data["Head"]["connections"][23]["weight"];
    params["RMDV_RMDD_ele"] = data["Head"]["connections"][24]["weight"];

    // Head neuron parameters.
    params["AIY_tau"] = data["Head"]["neurons"]["AIY"]["tau"];
    params["AIY_theta"] = data["Head"]["neurons"]["AIY"]["theta"];
    params["AWA_tau"] = data["Head"]["neurons"]["AWA"]["tau"];
    params["AWA_theta"] = data["Head"]["neurons"]["AWA"]["theta"];
    params["RIA_tau"] = data["Head"]["neurons"]["RIA"]["tau"];
    params["RIA_theta"] = data["Head"]["neurons"]["RIA"]["theta"];
    params["RMDD_tau"] = data["Head"]["neurons"]["RMDD"]["tau"];
    params["RMDD_theta"] = data["Head"]["neurons"]["RMDD"]["theta"];
    params["RMDV_tau"] = data["Head"]["neurons"]["RMDV"]["tau"];
    params["RMDV_theta"] = data["Head"]["neurons"]["RMDV"]["theta"];
    params["SMDD_tau"] = data["Head"]["neurons"]["SMDD"]["tau"];
    params["SMDD_theta"] = data["Head"]["neurons"]["SMDD"]["theta"];
    params["SMDV_tau"] = data["Head"]["neurons"]["SMDV"]["tau"];
    params["SMDV_theta"] = data["Head"]["neurons"]["SMDV"]["theta"];

    // NMJ parameters.
    params["NMJ_DB"] = data["NMJ"]["DB"];
    params["NMJ_DD"] = data["NMJ"]["DD"];
    params["NMJ_RMDD"] = data["NMJ"]["RMDD"];
    params["NMJ_RMDV"] = data["NMJ"]["RMDV"];
    params["NMJ_SMDD"] = data["NMJ"]["SMDD"];
    params["NMJ_SMDV"] = data["NMJ"]["SMDV"];
    params["NMJ_VBA"] = data["NMJ"]["VBA"];
    params["NMJ_VBP"] = data["NMJ"]["VBP"];
    params["NMJ_VDA"] = data["NMJ"]["VDA"];
    params["NMJ_VDP"] = data["NMJ"]["VDP"];

    // Stretch receptor parameters.
    params["sr_headgain"] = data["StretchReceptors"]["Head_gain"];
    params["sr_vcgain"] = data["StretchReceptors"]["VC_gain"];

    // VC connection parameters.
    params["DB_DB"] = data["VentralCord"]["connections"][0]["weight"];
    params["VBA_VBA"] = data["VentralCord"]["connections"][1]["weight"];
    params["VBP_VBP"] = data["VentralCord"]["connections"][2]["weight"];
    params["DD_DD"] = data["VentralCord"]["connections"][3]["weight"];
    params["VDA_VDA"] = data["VentralCord"]["connections"][4]["weight"];
    params["VDP_VDP"] = data["VentralCord"]["connections"][5]["weight"];
    params["DB_DD"] = data["VentralCord"]["connections"][6]["weight"];
    params["VBA_VDA"] = data["VentralCord"]["connections"][7]["weight"];
    params["VBP_VDP"] = data["VentralCord"]["connections"][8]["weight"];
    params["DB_VDA"] = data["VentralCord"]["connections"][9]["weight"];
    params["DB_VDP"] = data["VentralCord"]["connections"][10]["weight"];
    params["VBA_DD"] = data["VentralCord"]["connections"][11]["weight"];
    params["VBP_DD"] = data["VentralCord"]["connections"][12]["weight"];
    params["DD_VDA"] = data["VentralCord"]["connections"][13]["weight"];
    params["DD_VDA_ele"] = data["VentralCord"]["connections"][15]["weight"];
    params["DD_VDP_ele"] = data["VentralCord"]["connections"][16]["weight"];
    params["VDA_VDP_ele"] = data["VentralCord"]["connections"][17]["weight"];
    params["VBA_VBP_ele"] = data["VentralCord"]["connections"][18]["weight"];

    // Forward connection parameters.
    params["fwd_DB_DB"] = data["VentralCord"]["connections_fwd"][0]["weight"];
    params["fwd_VBP_VBA"] = data["VentralCord"]["connections_fwd"][1]["weight"];
    params["fwd_DD_DD"] = data["VentralCord"]["connections_fwd"][2]["weight"];
    params["fwd_VDP_VDA"] = data["VentralCord"]["connections_fwd"][3]["weight"];
    params["fwd_VBP_DB"] = data["VentralCord"]["connections_fwd"][4]["weight"];

    // VC neuron parameters.
    params["DB_tau"] = data["VentralCord"]["neurons"]["DB"]["tau"];
    params["DB_theta"] = data["VentralCord"]["neurons"]["DB"]["theta"];
    params["DD_tau"] = data["VentralCord"]["neurons"]["DD"]["tau"];
    params["DD_theta"] = data["VentralCord"]["neurons"]["DD"]["theta"];
    params["VBA_tau"] = data["VentralCord"]["neurons"]["VBA"]["tau"];
    params["VBA_theta"] = data["VentralCord"]["neurons"]["VBA"]["theta"];
    params["VBP_tau"] = data["VentralCord"]["neurons"]["VBP"]["tau"];
    params["VBP_theta"] = data["VentralCord"]["neurons"]["VDA"]["theta"];
    params["VDA_tau"] = data["VentralCord"]["neurons"]["VDA"]["tau"];
    params["VDA_theta"] = data["VentralCord"]["neurons"]["VDA"]["theta"];
    params["VDP_tau"] = data["VentralCord"]["neurons"]["VDP"]["tau"];
    params["VDP_theta"] = data["VentralCord"]["neurons"]["VDP"]["theta"];

    return params;
}

// Writes parameters to input/params.json.
int write_params(std::map<std::string, double>& params, const std::string& out) {

    alpha = -(params.find("alpha") -> second);; beta = params.find("beta") -> second;;
    gamma = params.find("gamma") -> second;; kappa = params.find("kappa") -> second;;
    lambda = params.find("lambda") -> second;;

    AWA_AIY = params.find("AWA_AIY") -> second;; AIY_AIY = params.find("AIY_AIY") -> second;;
    AIY_RIA = params.find("AIY_RIA") -> second;; RIA_RMDD = params.find("RIA_RMDD") -> second;;
    RIA_RMDV = params.find("RIA_RMDV") -> second;; SMDD_RIA = params.find("SMDD_RIA") -> second;;
    SMDV_RIA = params.find("SMDV_RIA") -> second;; RIA_RIA = params.find("RIA_RIA") -> second;;
    RIA_SMDD = params.find("RIA_SMDD") -> second;; RIA_SMDV = params.find("RIA_SMDV") -> second;;
    RMDD_RIA = params.find("RMDD_RIA") -> second;; RMDV_RIA = params.find("RMDV_RIA") -> second;;
    SMDD_SMDD = params.find("SMDD_SMDD") -> second;; SMDV_SMDV = params.find("SMDV_SMDV") -> second;;
    RMDD_RMDD = params.find("RMDD_RMDD") -> second;; RMDV_RMDV = params.find("RMDV_RMDV") -> second;;
    SMDD_SMDV = params.find("SMDD_SMDV") -> second;; SMDV_SMDD = params.find("SMDV_SMDD") -> second;;
    SMDD_RMDV = params.find("SMDD_RMDV") -> second;; SMDV_RMDD = params.find("SMDV_RMDD") -> second;;
    RMDD_RMDV = params.find("RMDD_RMDV") -> second;; RMDV_RMDD = params.find("RMDV_RMDD") -> second;;
    SMDD_RMDD_ele = params.find("SMDD_RMDD_ele") -> second;; SMDV_RMDV_ele = params.find("SMDV_RMDV_ele") -> second;;
    RMDV_RMDD_ele = params.find("RMDV_RMDD_ele") -> second;;

    AIY_tau = params.find("AIY_tau") -> second;; AIY_theta = params.find("AIY_theta") -> second;;
    AWA_tau = params.find("AWA_tau") -> second;; AWA_theta = params.find("AWA_theta") -> second;;
    RIA_tau = params.find("RIA_tau") -> second;; RIA_theta = params.find("RIA_theta") -> second;;
    RMDD_tau = params.find("RMDD_tau") -> second;; RMDD_theta = params.find("RMDD_theta") -> second;;
    RMDV_tau = params.find("RMDV_tau") -> second;; RMDV_theta = params.find("RMDV_theta") -> second;;
    SMDD_tau = params.find("SMDD_tau") -> second;; SMDD_theta = params.find("SMDD_theta") -> second;;
    SMDV_tau = params.find("SMDV_tau") -> second;; SMDV_theta = params.find("SMDV_theta") -> second;;

    NMJ_DB = params.find("NMJ_DB") -> second;; NMJ_DD = params.find("NMJ_DD") -> second;;
    NMJ_RMDD = params.find("NMJ_RMDD") -> second;; NMJ_RMDV = params.find("NMJ_RMDV") -> second;;
    NMJ_SMDD =params.find("NMJ_SMDD") -> second;; NMJ_SMDV = params.find("NMJ_SMDV") -> second;;
    NMJ_VBA = params.find("NMJ_VBA") -> second;; NMJ_VBP = params.find("NMJ_VBP") -> second;;
    NMJ_VDA = params.find("NMJ_VDA") -> second;; NMJ_VDP = params.find("NMJ_VDP") -> second;;

    sr_headgain = params.find("sr_headgain") -> second;; sr_vcgain = params.find("sr_vcgain") -> second;;

    DB_DB = params.find("DB_DB") -> second; VBA_VBA = params.find("VBA_VBA") -> second;
    VBP_VBP = params.find("VBP_VBP") -> second;; DD_DD = params.find("DD_DD") -> second;
    VDA_VDA = params.find("VDA_VDA") -> second; VDP_VDP = params.find("VDP_VDP") -> second;
    DB_DD = params.find("DB_DD") -> second; VBA_VDA = params.find("VBA_VDA") -> second;
    VBP_VDP = params.find("VBP_VDP") -> second; DB_VDA = params.find("DB_VDA") -> second;
    DB_VDP = params.find("DB_VDP") -> second; VBA_DD = params.find("VBA_DD") -> second;
    VBP_DD = params.find("VBP_DD") -> second; DD_VDA = params.find("DD_VDA") -> second;
    DD_VDA_ele = params.find("DD_VDA_ele") -> second; DD_VDP_ele = params.find("DD_VDP_ele") -> second;
    VDA_VDP_ele = params.find("VDA_VDP_ele") -> second; VBA_VBP_ele = params.find("VBA_VBP_ele") -> second;

    fwd_DB_DB = params.find("fwd_DB_DB") -> second; fwd_VBP_VBA = params.find("fwd_VBP_VBA") -> second;
    fwd_DD_DD = params.find("fwd_DD_DD") -> second; fwd_VDP_VDA = params.find("fwd_VDP_VDA") -> second;
    fwd_VBP_DB = params.find("fwd_VBP_DB") -> second;

    DB_tau = params.find("DB_tau") -> second; DB_theta = params.find("DB_theta") -> second;
    DD_tau = params.find("DD_tau") -> second; DD_theta = params.find("DD_theta") -> second;
    VBA_tau = params.find("VBA_tau") -> second; VBA_theta = params.find("VBA_theta") -> second;
    VBP_tau = params.find("VBP_tau") -> second; VBP_theta = params.find("VBP_theta") -> second;
    VDA_tau = params.find("VDA_tau") -> second; VDA_theta = params.find("VDA_theta") -> second;
    VDP_tau = params.find("VDP_tau") -> second;; VDP_theta = params.find("VDP_theta") -> second;

    std::ifstream f("input/params.json");
    json data;
    f >> data;

    data["ChemoReceptors"]["alpha"] = alpha;
    data["ChemoReceptors"]["beta"] = beta;
    data["ChemoReceptors"]["gamma"] = p_gamma;
    data["ChemoReceptors"]["kappa"] = kappa;
    data["ChemoReceptors"]["lambda"] = lambda;

    data["Head"]["connections"][0]["weight"] = AWA_AIY;
    data["Head"]["connections"][1]["weight"] = AIY_AIY;
    data["Head"]["connections"][2]["weight"] = AIY_RIA;
    data["Head"]["connections"][3]["weight"] = RIA_RMDD;
    data["Head"]["connections"][4]["weight"] = RIA_RMDV;
    data["Head"]["connections"][5]["weight"] = SMDD_RIA;
    data["Head"]["connections"][6]["weight"] = SMDV_RIA;
    data["Head"]["connections"][7]["weight"] = RIA_RIA;
    data["Head"]["connections"][8]["weight"] = RIA_SMDD;
    data["Head"]["connections"][9]["weight"] = RIA_SMDV;
    data["Head"]["connections"][10]["weight"] = RMDD_RIA;
    data["Head"]["connections"][11]["weight"] = RMDV_RIA;
    data["Head"]["connections"][12]["weight"] = SMDD_SMDD;
    data["Head"]["connections"][13]["weight"] = SMDV_SMDV;
    data["Head"]["connections"][14]["weight"] = RMDD_RMDD;
    data["Head"]["connections"][15]["weight"] = RMDV_RMDV;
    data["Head"]["connections"][16]["weight"] = SMDD_SMDV;
    data["Head"]["connections"][17]["weight"] = SMDV_SMDD;
    data["Head"]["connections"][18]["weight"] = SMDD_RMDV;
    data["Head"]["connections"][19]["weight"] = SMDV_RMDD;
    data["Head"]["connections"][20]["weight"] = RMDD_RMDV;
    data["Head"]["connections"][21]["weight"] = RMDV_RMDD;
    data["Head"]["connections"][22]["weight"] = SMDD_RMDD_ele;
    data["Head"]["connections"][23]["weight"] = SMDV_RMDV_ele;
    data["Head"]["connections"][24]["weight"] = RMDV_RMDD_ele;

    data["Head"]["neurons"]["AIY"]["tau"] = AIY_tau;
    data["Head"]["neurons"]["AIY"]["theta"] = AIY_theta;
    data["Head"]["neurons"]["AWA"]["tau"] = AWA_tau;
    data["Head"]["neurons"]["AWA"]["theta"] = AWA_theta;
    data["Head"]["neurons"]["RIA"]["tau"] = RIA_tau;
    data["Head"]["neurons"]["RIA"]["theta"] = RIA_theta;
    data["Head"]["neurons"]["RMDD"]["tau"] = RMDD_tau;
    data["Head"]["neurons"]["RMDD"]["theta"] = RMDD_theta;
    data["Head"]["neurons"]["RMDV"]["tau"] = RMDV_tau;
    data["Head"]["neurons"]["RMDV"]["theta"] = RMDV_theta;
    data["Head"]["neurons"]["SMDD"]["tau"] = SMDD_tau;
    data["Head"]["neurons"]["SMDD"]["theta"] = SMDD_theta;
    data["Head"]["neurons"]["SMDV"]["tau"] = SMDV_tau;
    data["Head"]["neurons"]["SMDV"]["theta"] = SMDV_theta;

    data["NMJ"]["DB"] = NMJ_DB;
    data["NMJ"]["DD"] = NMJ_DD;
    data["NMJ"]["RMDD"] = NMJ_RMDD;
    data["NMJ"]["RMDV"] = NMJ_RMDV;
    data["NMJ"]["SMDD"] = NMJ_SMDD;
    data["NMJ"]["SMDV"] = NMJ_SMDV;
    data["NMJ"]["VBA"] = NMJ_VBA;
    data["NMJ"]["VBP"] = NMJ_VBP;
    data["NMJ"]["VDA"] = NMJ_VDA;
    data["NMJ"]["VDP"] = NMJ_VDP;

    data["StretchReceptors"]["Head_gain"] = sr_headgain;
    data["StretchReceptors"]["VC_gain"] = sr_vcgain;

    data["VentralCord"]["connections"][0]["weight"] = DB_DB;
    data["VentralCord"]["connections"][1]["weight"] = VBA_VBA;
    data["VentralCord"]["connections"][2]["weight"] = VBP_VBP;
    data["VentralCord"]["connections"][3]["weight"] = DD_DD;
    data["VentralCord"]["connections"][4]["weight"] = VDA_VDA;
    data["VentralCord"]["connections"][5]["weight"] = VDP_VDP;
    data["VentralCord"]["connections"][6]["weight"] = DB_DD;
    data["VentralCord"]["connections"][7]["weight"] = VBA_VDA;
    data["VentralCord"]["connections"][8]["weight"] = VBP_VDP;
    data["VentralCord"]["connections"][9]["weight"] = DB_VDA;
    data["VentralCord"]["connections"][10]["weight"] = DB_VDP;
    data["VentralCord"]["connections"][11]["weight"] = VBA_DD;
    data["VentralCord"]["connections"][12]["weight"] = VBP_DD;
    data["VentralCord"]["connections"][13]["weight"] = DD_VDA;
    data["VentralCord"]["connections"][14]["weight"] = DD_VDA;
    data["VentralCord"]["connections"][15]["weight"] = DD_VDA_ele;
    data["VentralCord"]["connections"][16]["weight"] = DD_VDP_ele;
    data["VentralCord"]["connections"][17]["weight"] = VDA_VDP_ele;
    data["VentralCord"]["connections"][18]["weight"] = VBA_VBP_ele;

    data["VentralCord"]["connections_fwd"][0]["weight"] = fwd_DB_DB;
    data["VentralCord"]["connections_fwd"][1]["weight"] = fwd_VBP_VBA;
    data["VentralCord"]["connections_fwd"][2]["weight"] = fwd_DD_DD;
    data["VentralCord"]["connections_fwd"][3]["weight"] = fwd_VDP_VDA;
    data["VentralCord"]["connections_fwd"][4]["weight"] = fwd_VBP_DB;

    data["VentralCord"]["neurons"]["DB"]["tau"] = DB_tau;
    data["VentralCord"]["neurons"]["DB"]["theta"] = DB_theta;
    data["VentralCord"]["neurons"]["DD"]["tau"] = DD_tau;
    data["VentralCord"]["neurons"]["DD"]["theta"] = DD_theta;
    data["VentralCord"]["neurons"]["VBA"]["tau"] = VBA_tau;
    data["VentralCord"]["neurons"]["VBA"]["theta"] = VBA_theta;
    data["VentralCord"]["neurons"]["VBP"]["tau"] = VBP_tau;
    data["VentralCord"]["neurons"]["VBP"]["theta"] = VBP_theta;
    data["VentralCord"]["neurons"]["VDA"]["theta"] = VBP_theta;
    data["VentralCord"]["neurons"]["VDA"]["tau"] = VDA_tau;
    data["VentralCord"]["neurons"]["VDA"]["theta"] = VDA_theta;
    data["VentralCord"]["neurons"]["VDP"]["tau"] = VDP_tau;
    data["VentralCord"]["neurons"]["VDP"]["theta"] = VDP_theta;

    data["simulation"]["duration"] = duration;
    data["simulation"]["angle"] = angle;

    std::ofstream g(out);

    if (g.is_open()) {
        g << data.dump(4);
        std::cout << "JSON data has been written to " << out << std::endl;
    } else {
        std::cerr << "Failed to open the file for writing." << std::endl;
    }
    g.close();
    return 0;
};

// Write base params (for hybrid initialization). 
int write_base_params() {
    alpha = 4.0; beta = 15.0;
    foodPos_x = 0.002; foodPos_y = 0.003;
    gamma = 2.0; kappa = 244.00503969854154;
    lambda = -54032.94745512404;

    AWA_AIY = -12.300635881587867; AIY_AIY = 0.26166624096194646;
    AIY_RIA = -45.91277074924063; RIA_RMDD = -16.07013704852836;
    RIA_RMDV = -16.07013704852836; SMDD_RIA = 2.6014517231130534;
    SMDV_RIA = 2.6014517231130534; RIA_RIA = -0.06285081961384109;
    RIA_SMDD = 30.450189759537384; RIA_SMDV = 30.450189759537384;
    RMDD_RIA = 0.7417008490490758; RMDV_RIA = 0.7417008490490758;
    SMDD_SMDD = -14.9121; SMDV_SMDV = -14.9121;
    RMDD_RMDD = 6.62512; RMDV_RMDV = 6.62512;
    SMDD_SMDV = -11.2755; SMDV_SMDD = -11.2755;
    SMDD_RMDV = 14.9933; SMDV_RMDD = 14.9933;
    RMDD_RMDV = -11.6075; RMDV_RMDD = -11.6075;
    SMDD_RMDD_ele = 0.0199558; SMDV_RMDV_ele = 0.0199558;
    RMDV_RMDD_ele = 1.51095;

    AIY_tau = 0.01; AIY_theta = -0.7247370228031592;
    AWA_tau = 0.005; AWA_theta = -3.2532581458769427;
    RIA_tau = 0.01; RIA_theta = -2.1267449437585038;
    RMDD_tau = 0.502018; RMDD_theta = -3.65014;
    RMDV_tau = 0.502018; RMDV_theta = -3.65014;
    SMDD_tau = 0.500646; SMDD_theta = 6.27095;
    SMDV_tau = 0.500646; SMDV_theta = 6.27095;

    NMJ_DB = 1.0; NMJ_DD = -0.00105857;
    NMJ_RMDD = 1.0; NMJ_RMDV = 1.0;
    NMJ_SMDD = 0.000371977; NMJ_SMDV = 0.000371977;
    NMJ_VBA = 1.0; NMJ_VBP = 1.0;
    NMJ_VDA = -0.00105857; NMJ_VDP = -0.00105857;

    sr_headgain = -120.535; sr_vcgain = -165.138;

    DB_DB = -14.9544; VBA_VBA = -14.9544;
    VBP_VBP = -14.9544; DD_DD = 1.13632;
    VDA_VDA = 1.13632; VDP_VDP = 1.13632;
    DB_DD = -9.79669; VBA_VDA = -9.79669;
    VBP_VDP = -9.79669; DB_VDA = 12.3531;
    DB_VDP = 12.3531; VBA_DD = 6.17657;
    VBP_DD = 6.17657; DD_VDA = 1.66374;
    DD_VDA_ele = 1.46108; DD_VDP_ele = 1.46108;
    VDA_VDP_ele = 0.994941; VBA_VBP_ele = 0.569164;

    fwd_DB_DB = 0.569164; fwd_VBP_VBA = 0.569164;
    fwd_DD_DD = 0.994941; fwd_VDP_VDA = 0.994941;
    fwd_VBP_DB = 1.10313;

    DB_tau = 0.502865; DB_theta = 4.80699;
    DD_tau = 1.67861; DD_theta = -2.94785;
    VBA_tau = 0.502865; VBA_theta = 4.80699;
    VBP_tau = 0.502865; VBP_theta = 4.80699;
    VDA_tau = 1.67861; VDA_theta = -2.94785;
    VDP_tau = 1.67861; VDP_theta = -2.94785;

    std::ifstream f("input/params.json");
    json data;
    f >> data;

    data["ChemoReceptors"]["alpha"] = alpha;
    data["ChemoReceptors"]["beta"] = beta;
    data["ChemoReceptors"]["foodPos"]["x"] = foodPos_x;
    data["ChemoReceptors"]["foodPos"]["y"] = foodPos_y;
    data["ChemoReceptors"]["gamma"] = gamma;
    data["ChemoReceptors"]["kappa"] = kappa;
    data["ChemoReceptors"]["lambda"] = lambda;

    data["Head"]["connections"][0]["weight"] = AWA_AIY;
    data["Head"]["connections"][1]["weight"] = AIY_AIY;
    data["Head"]["connections"][2]["weight"] = AIY_RIA;
    data["Head"]["connections"][3]["weight"] = RIA_RMDD;
    data["Head"]["connections"][4]["weight"] = RIA_RMDV;
    data["Head"]["connections"][5]["weight"] = SMDD_RIA;
    data["Head"]["connections"][6]["weight"] = SMDV_RIA;
    data["Head"]["connections"][7]["weight"] = RIA_RIA;
    data["Head"]["connections"][8]["weight"] = RIA_SMDD;
    data["Head"]["connections"][9]["weight"] = RIA_SMDV;
    data["Head"]["connections"][10]["weight"] = RMDD_RIA;
    data["Head"]["connections"][11]["weight"] = RMDV_RIA;
    data["Head"]["connections"][12]["weight"] = SMDD_SMDD;
    data["Head"]["connections"][13]["weight"] = SMDV_SMDV;
    data["Head"]["connections"][14]["weight"] = RMDD_RMDD;
    data["Head"]["connections"][15]["weight"] = RMDV_RMDV;
    data["Head"]["connections"][16]["weight"] = SMDD_SMDV;
    data["Head"]["connections"][17]["weight"] = SMDV_SMDD;
    data["Head"]["connections"][18]["weight"] = SMDD_RMDV;
    data["Head"]["connections"][19]["weight"] = SMDV_RMDD;
    data["Head"]["connections"][20]["weight"] = RMDD_RMDV;
    data["Head"]["connections"][21]["weight"] = RMDV_RMDD;
    data["Head"]["connections"][22]["weight"] = SMDD_RMDD_ele;
    data["Head"]["connections"][23]["weight"] = SMDV_RMDV_ele;
    data["Head"]["connections"][24]["weight"] = RMDV_RMDD_ele;

    data["Head"]["neurons"]["AIY"]["tau"] = AIY_tau;
    data["Head"]["neurons"]["AIY"]["theta"] = AIY_theta;
    data["Head"]["neurons"]["AWA"]["tau"] = AWA_tau;
    data["Head"]["neurons"]["AWA"]["theta"] = AWA_theta;
    data["Head"]["neurons"]["RIA"]["tau"] = RIA_tau;
    data["Head"]["neurons"]["RIA"]["theta"] = RIA_theta;
    data["Head"]["neurons"]["RMDD"]["tau"] = RMDD_tau;
    data["Head"]["neurons"]["RMDD"]["theta"] = RMDD_theta;
    data["Head"]["neurons"]["RMDV"]["tau"] = RMDV_tau;
    data["Head"]["neurons"]["RMDV"]["theta"] = RMDV_theta;
    data["Head"]["neurons"]["SMDD"]["tau"] = SMDD_tau;
    data["Head"]["neurons"]["SMDD"]["theta"] = SMDD_theta;
    data["Head"]["neurons"]["SMDV"]["tau"] = SMDV_tau;
    data["Head"]["neurons"]["SMDV"]["theta"] = SMDV_theta;

    data["NMJ"]["DB"] = NMJ_DB;
    data["NMJ"]["DD"] = NMJ_DD;
    data["NMJ"]["RMDD"] = NMJ_RMDD;
    data["NMJ"]["RMDV"] = NMJ_RMDV;
    data["NMJ"]["SMDD"] = NMJ_SMDD;
    data["NMJ"]["SMDV"] = NMJ_SMDV;
    data["NMJ"]["VBA"] = NMJ_VBA;
    data["NMJ"]["VBP"] = NMJ_VBP;
    data["NMJ"]["VDA"] = NMJ_VDA;
    data["NMJ"]["VDP"] = NMJ_VDP;

    data["StretchReceptors"]["Head_gain"] = sr_headgain;
    data["StretchReceptors"]["VC_gain"] = sr_vcgain;

    data["VentralCord"]["connections"][0]["weight"] = DB_DB;
    data["VentralCord"]["connections"][1]["weight"] = VBA_VBA;
    data["VentralCord"]["connections"][2]["weight"] = VBP_VBP;
    data["VentralCord"]["connections"][3]["weight"] = DD_DD;
    data["VentralCord"]["connections"][4]["weight"] = VDA_VDA;
    data["VentralCord"]["connections"][5]["weight"] = VDP_VDP;
    data["VentralCord"]["connections"][6]["weight"] = DB_DD;
    data["VentralCord"]["connections"][7]["weight"] = VBA_VDA;
    data["VentralCord"]["connections"][8]["weight"] = VBP_VDP;
    data["VentralCord"]["connections"][9]["weight"] = DB_VDA;
    data["VentralCord"]["connections"][10]["weight"] = DB_VDP;
    data["VentralCord"]["connections"][11]["weight"] = VBA_DD;
    data["VentralCord"]["connections"][12]["weight"] = VBP_DD;
    data["VentralCord"]["connections"][13]["weight"] = DD_VDA;
    data["VentralCord"]["connections"][15]["weight"] = DD_VDA_ele;
    data["VentralCord"]["connections"][16]["weight"] = DD_VDP_ele;
    data["VentralCord"]["connections"][17]["weight"] = VDA_VDP_ele;
    data["VentralCord"]["connections"][18]["weight"] = VBA_VBP_ele;

    data["VentralCord"]["connections_fwd"][0]["weight"] = fwd_DB_DB;
    data["VentralCord"]["connections_fwd"][1]["weight"] = fwd_VBP_VBA;
    data["VentralCord"]["connections_fwd"][2]["weight"] = fwd_DD_DD;
    data["VentralCord"]["connections_fwd"][3]["weight"] = fwd_VDP_VDA;
    data["VentralCord"]["connections_fwd"][4]["weight"] = fwd_VBP_DB;

    data["VentralCord"]["neurons"]["DB"]["tau"] = DB_tau;
    data["VentralCord"]["neurons"]["DB"]["theta"] = DB_theta;
    data["VentralCord"]["neurons"]["DD"]["tau"] = DD_tau;
    data["VentralCord"]["neurons"]["DD"]["theta"] = DD_theta;
    data["VentralCord"]["neurons"]["VBA"]["tau"] = VBA_tau;
    data["VentralCord"]["neurons"]["VBA"]["theta"] = VBA_theta;
    data["VentralCord"]["neurons"]["VBP"]["tau"] = VBP_tau;
    data["VentralCord"]["neurons"]["VDA"]["theta"] = VBP_theta;
    data["VentralCord"]["neurons"]["VDA"]["tau"] = VDA_tau;
    data["VentralCord"]["neurons"]["VDA"]["theta"] = VDA_theta;
    data["VentralCord"]["neurons"]["VDP"]["tau"] = VDP_tau;
    data["VentralCord"]["neurons"]["VDP"]["theta"] = VDP_theta;

    std::ofstream base_params_file("input/params.json");

    if (base_params_file.is_open()) {
        base_params_file << data.dump(4);
        std::cout << "JSON data has been written to input/params.json" << std::endl;
    } else {
        std::cerr << "Failed to open the file for writing." << std::endl;
    }
    base_params_file.close();
    return 0;
}

// Write random parameters.
json write_random_params() {

    std::srand(std::time(0));

 		alpha = gen.get_rand(0.4, 40); beta = gen.get_rand(5, 150);
    	foodPos_x = 0.002; foodPos_y = 0.003;
    	p_gamma = gen.get_rand(0.2, 20);

    	AWA_AIY = gen.get_rand(-120, -1.2); AIY_AIY = gen.get_rand(0.02, 2);
    	AIY_RIA = gen.get_rand(-450, -4.5);
   	 	double RIA_RMD_rand = gen.get_rand(-160, -1.6);
    	RIA_RMDD = RIA_RMD_rand; RIA_RMDV = RIA_RMD_rand;
    	double SMD_RIA_rand = gen.get_rand(0.2, 20);
    	SMDD_RIA = SMD_RIA_rand; SMDV_RIA = SMD_RIA_rand;
    	RIA_RIA = gen.get_rand(-0.006, -0.6);
    	double RIA_SMD_rand = gen.get_rand(3, 300);
    	RIA_SMDD = RIA_SMD_rand; RIA_SMDV = RIA_SMD_rand;
    	double RMD_RIA_rand = gen.get_rand(0.07, 7);
    	RMDD_RIA = RMD_RIA_rand; RMDV_RIA = RMD_RIA_rand;
    	double a_SMD_SMD_rand = gen.get_rand(-28.0, -7.0);
    	SMDD_SMDD = a_SMD_SMD_rand; SMDV_SMDV = a_SMD_SMD_rand;
    	double RMD_RMD_rand = gen.get_rand(3, 12);
    	RMDD_RMDD = RMD_RMD_rand; RMDV_RMDV = RMD_RMD_rand;
    	double b_SMD_SMD_rand = gen.get_rand(-22, -5);
    	SMDD_SMDV = b_SMD_SMD_rand; SMDV_SMDD = b_SMD_SMD_rand;
    	double SMD_RMD_rand = gen.get_rand(7, 28);
    	SMDD_RMDV = SMD_RMD_rand; SMDV_RMDD = SMD_RMD_rand;
    	double b_RMD_RMD_rand = gen.get_rand(-22, -5);
    	RMDD_RMDV = b_RMD_RMD_rand; RMDV_RMDD = b_RMD_RMD_rand;
    	double SMD_RMD_ele_rand = gen.get_rand(0.005, 0.02);
    	SMDD_RMDD_ele = std::abs(SMD_RMD_ele_rand); SMDV_RMDV_ele = std::abs(SMD_RMD_ele_rand);
    	RMDV_RMDD_ele = gen.get_rand(0.7, 3);

    	AIY_tau = gen.get_rand(0.01, 1); AIY_theta = gen.get_rand(-7, -0.07);
    	AWA_tau = gen.get_rand(0.01, 1); AWA_theta = gen.get_rand(-30, -0.3);
    	RIA_tau = gen.get_rand(0.01, 1); RIA_theta = gen.get_rand(-20, -0.2);
    	double RMD_tau_rand = gen.get_rand(0.1, 1);
    	RMDD_tau = RMD_tau_rand; RMDV_tau = RMD_tau_rand;
    	double RMD_theta_rand = gen.get_rand(-6,-1.5);
    	RMDD_theta = RMD_theta_rand; RMDV_theta = RMD_theta_rand;
    	double SMD_tau_rand = gen.get_rand(0.1, 1);
    	SMDD_tau = SMD_tau_rand; SMDV_tau = SMD_tau_rand;
    	double SMD_theta_rand = gen.get_rand(3, 12);
    	SMDV_theta = SMD_theta_rand; SMDV_theta = SMD_theta_rand;

        std::ifstream f("input/params.json");
        json data;
        f >> data;
    
        data["ChemoReceptors"]["alpha"] = alpha;
        data["ChemoReceptors"]["beta"] = beta;
        data["ChemoReceptors"]["foodPos"]["x"] = foodPos_x;
        data["ChemoReceptors"]["gamma"] = p_gamma;
        data["ChemoReceptors"]["kappa"] = kappa;
        data["ChemoReceptors"]["lambda"] = lambda;
        data["ChemoReceptors"]["foodPos"]["x"] = foodPos_x;
        data["ChemoReceptors"]["foodPos"]["y"] = foodPos_y;
    
        data["Head"]["connections"][0]["weight"] = AWA_AIY;
        data["Head"]["connections"][1]["weight"] = AIY_AIY;
        data["Head"]["connections"][2]["weight"] = AIY_RIA;
        data["Head"]["connections"][3]["weight"] = RIA_RMDD;
        data["Head"]["connections"][4]["weight"] = RIA_RMDV;
        data["Head"]["connections"][5]["weight"] = SMDD_RIA;
        data["Head"]["connections"][6]["weight"] = SMDV_RIA;
        data["Head"]["connections"][7]["weight"] = RIA_RIA;
        data["Head"]["connections"][8]["weight"] = RIA_SMDD;
        data["Head"]["connections"][9]["weight"] = RIA_SMDV;
        data["Head"]["connections"][10]["weight"] = RMDD_RIA;
        data["Head"]["connections"][11]["weight"] = RMDV_RIA;
        data["Head"]["connections"][12]["weight"] = SMDD_SMDD;
        data["Head"]["connections"][13]["weight"] = SMDV_SMDV;
        data["Head"]["connections"][14]["weight"] = RMDD_RMDD;
        data["Head"]["connections"][15]["weight"] = RMDV_RMDV;
        data["Head"]["connections"][16]["weight"] = SMDD_SMDV;
        data["Head"]["connections"][17]["weight"] = SMDV_SMDD;
        data["Head"]["connections"][18]["weight"] = SMDD_RMDV;
        data["Head"]["connections"][19]["weight"] = SMDV_RMDD;
        data["Head"]["connections"][20]["weight"] = RMDD_RMDV;
        data["Head"]["connections"][21]["weight"] = RMDV_RMDD;
        data["Head"]["connections"][22]["weight"] = SMDD_RMDD_ele;
        data["Head"]["connections"][23]["weight"] = SMDV_RMDV_ele;
        data["Head"]["connections"][24]["weight"] = RMDV_RMDD_ele;
    
        data["Head"]["neurons"]["AIY"]["tau"] = AIY_tau;
        data["Head"]["neurons"]["AIY"]["theta"] = AIY_theta;
        data["Head"]["neurons"]["AWA"]["tau"] = AWA_tau;
        data["Head"]["neurons"]["AWA"]["theta"] = AWA_theta;
        data["Head"]["neurons"]["RIA"]["tau"] = RIA_tau;
        data["Head"]["neurons"]["RIA"]["theta"] = RIA_theta;
        data["Head"]["neurons"]["RMDD"]["tau"] = RMDD_tau;
        data["Head"]["neurons"]["RMDD"]["theta"] = RMDD_theta;
        data["Head"]["neurons"]["RMDV"]["tau"] = RMDV_tau;
        data["Head"]["neurons"]["RMDV"]["theta"] = RMDV_theta;
        data["Head"]["neurons"]["SMDD"]["tau"] = SMDD_tau;
        data["Head"]["neurons"]["SMDD"]["theta"] = SMDD_theta;
        data["Head"]["neurons"]["SMDV"]["tau"] = SMDV_tau;
        data["Head"]["neurons"]["SMDV"]["theta"] = SMDV_theta;
    
        std::ofstream random_params_file("input/params.json");
        if (random_params_file.is_open())
        {
            random_params_file << data.dump(4);
            std::cout << "JSON data has been written to input/params.json" << std::endl;
        } else
        {
            std::cerr << "Failed to open the file for writing." << std::endl;
        }
        random_params_file.close();
    
        return data;
    }
    return data;
};}
