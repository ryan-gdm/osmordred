#include <RDBoost/python.h>
#include <Osmordred/Osmordred.h>

namespace python = boost::python;

using namespace RDKit;

BOOST_PYTHON_MODULE(rdOsmordred) {
    python::scope().attr("__doc__") = "Functions for Mordred calculations";

    python::def("CalcABCIndex", calcABCIndex,
        "CalcABCIndex Lorem Ipsum\n");
    python::def("CalcAcidBase", calcAcidBase,
        "CalcAcidBase Lorem Ipsum\n");
    python::def("CalcAromatic", calcAromatic,
        "CalcAromatic Lorem Ipsum\n");
    python::def("CalcAtomCount", calcAtomCounts,
        "CalcAtomCounts Lorem Ipsum\n");
    python::def("CalcBalabanJ", calcBalabanJ,
        "CalcBalabanJ Lorem Ipsum\n");
    python::def("CalcBertzCT", calcBertzCT,
        "CalcBertzCT Lorem Ipsum\n");
    python::def("CalcBondCount", calcBondCounts,
        "CalcBondCounts Lorem Ipsum\n");
    python::def("CalcVertexAdjacencyInformation", calcVertexAdjacencyInformation,
        "CalcVertexAdjacencyInformation Lorem Ipsum\n");
    python::def("CalcWeight", calcWeight,
        "CalcWeight Lorem Ipsum\n");
    python::def("CalcWienerIndex", calcWienerIndex,
        "CalcWienerIndex Lorem Ipsum\n");
    python::def("CalcVdwVolumeABC", calcVdwVolumeABC,
        "CalcVdwVolumeABC Lorem Ipsum\n");
    python::def("CalcTopoPSA", calcTopoPSA,
        "CalcTopoPSA Lorem Ipsum\n");
    python::def("CalcSLogP", calcSLogP,
        "CalcSLogP Lorem Ipsum\n");
    python::def("CalcHydrogenBond", calcHydrogenBond,
        "CalcHydrogenBond Lorem Ipsum\n");
    python::def("CalcLogS", calcLogS,
        "CalcLogS Lorem Ipsum\n");
    python::def("CalcLipinski", calcLipinskiGhose,
        "CalcLipinskiGhose Lorem Ipsum\n");
    python::def("CalcMcGowanVolume", calcMcGowanVolume,
        "CalcMcGowanVolume Lorem Ipsum\n");
    python::def("CalcPolarizability", calcPolarizability,
        "CalcPolarizability Lorem Ipsum\n");
    python::def("CalcRotatableBond", calcRotatableBond,
        "CalcRotatableBond Lorem Ipsum\n");
    python::def("CalcFragmentComplexity", calcFragmentComplexity,
        "CalcFragmentComplexity Lorem Ipsum\n");
    python::def("CalcConstitutional", calcConstitutional,
        "CalcConstitutional Lorem Ipsum\n");
    python::def("CalcTopologicalIndex", calcTopologicalIndex,
        "CalcTopologicalIndex Lorem Ipsum\n");
    python::def("CalcDetourMatrixEigen", calcDetourMatrixDescs,
        "CalcDetourMatrixDescs Lorem Ipsum\n");
    python::def("CalcDetourMatrix", calcDetourMatrixDescsL,
        "CalcDetourMatrixDescsL Lorem Ipsum\n");
    python::def("CalcDistanceMatrixEigen", calcDistMatrixDescs,
        "CalcDistMatrixDescs Lorem Ipsum\n");
    python::def("CalcDistanceMatrix", calcDistMatrixDescsL,
        "CalcDistMatrixDescsL Lorem Ipsum\n");
    python::def("CalcAdjacencyMatrixEigen", calcAdjMatrixDescs,
        "CalcAdjMatrixDescs Lorem Ipsum\n");
    python::def("CalcAdjacencyMatrix", calcAdjMatrixDescsL,
        "CalcAdjMatrixDescsL Lorem Ipsum\n");
    python::def("CalcCarbonTypes", calcCarbonTypes,
        "CalcCarbonTypes Lorem Ipsum\n");
    python::def("CalcEccentricConnectivityIndex", calcEccentricConnectivityIndex,
        "CalcEccentricConnectivityIndex Lorem Ipsum\n");
    python::def("CalcBaryszMatrix", calcBaryszMatrixDescsL,
        "CalcBaryszMatrixDescsL Lorem Ipsum\n");
    python::def("CalcBaryszMatrixEigen", calcBaryszMatrixDescs,
        "CalcBaryszMatrixDescs Lorem Ipsum\n");
    python::def("CalcZagrebIndex", calcZagrebIndex,
        "CalcZagrebIndex Lorem Ipsum\n");
    python::def("CalcMoeType", calcMoeType,
        "CalcMoeType Lorem Ipsum\n");
    python::def("CalcMolecularDistanceEdge", calcMolecularDistanceEdgeDescs,
        "CalcMolecularDistanceEdgeDescs Lorem Ipsum\n");
    python::def("CalcEState", calcEStateDescs,
        "CalcEStateDescs Lorem Ipsum\n");
    python::def("CalcWalkCount", calcWalkCounts,
        "CalcWalkCounts Lorem Ipsum\n");
    python::def("CalcTopologicalCharge", calcTopologicalChargeDescs,
        "CalcTopologicalChargeDescs Lorem Ipsum\n");
    python::def("CalcChi", calcAllChiDescriptors,
        "CalcAllChiDescriptors Lorem Ipsum\n");
    python::def("CalcPathCount", calcPathCount,
        "CalcPathCount Lorem Ipsum\n");
    python::def("CalcKappaShapeIndex", calcKappaShapeIndex,
        "CalcKappaShapeIndex Lorem Ipsum\n");
    python::def("CalcRingCount", calcRingCount,
        "CalcRingCount Lorem Ipsum\n");
    python::def("CalcMolecularId", calcMolecularId,
        "CalcMolecularId Lorem Ipsum\n");
    python::def("CalcBCUT", calcBCUTs,
        "CalcBCUTs Lorem Ipsum\n");
    python::def("CalcAutocorrelation", calcAutoCorrelation,
        "CalcAutoCorrelation Lorem Ipsum\n");
    python::def("CalcFramework", calcFramework,
        "CalcFramework Lorem Ipsum\n");
    python::def("CalcExtendedTopochemicalAtom", calcExtendedTopochemicalAtom,
        "CalcExtendedTopochemicalAtom Lorem Ipsum\n");
    python::def("CalcExtendedTopochemicalAtom2", calculateETADescriptors,
        "CalculateETADescriptors Lorem Ipsum\n");
    python::def("CalcChipath", calcChipath,
        "CalcChipath Lorem Ipsum\n");
    python::def("CalcChichain", calcChichain,
        "CalcChichain Lorem Ipsum\n");
    python::def("CalcChicluster", calcChicluster,
        "CalcChicluster Lorem Ipsum\n");
    python::def("CalcChipathcluster", calcChipathcluster,
        "CalcChipathcluster Lorem Ipsum\n");
    python::def("CalcAcidicGroupCount", calcAcidicGroupCount,
        "CalcAcidicGroupCount Lorem Ipsum\n");
    python::def("CalcBasicGroupCount", calcBasicGroupCount,
        "CalcBasicGroupCount Lorem Ipsum\n");
    python::def("CalcCountAromaticAtoms", countAromaticAtoms,
        "CalcCountAromaticAtoms Lorem Ipsum");
    python::def("CalcCountAromaticBonds", countAromaticBonds,
        "CalcCountAromaticBonds Lorem Ipsum");
    python::def("CalcBEState", calcBEStateDescs,
        "CalcBEStateDescs Lorem Ipsum\n");
    python::def("CalcHEState", calcHEStateDescs,
        "CalcHEStateDescs Lorem Ipsum\n");
    python::def("CalcAlphaKappaShapeIndex", calcAlphaKappaShapeIndex,
        "CalcAlphaKappaShapeIndex Lorem Ipsum\n");
    python::def("CalcAbrahams", calcAbrahams,
        "CalcAbrahams Lorem Ipsum\n");
    python::def("CalcPol", calcPol,
        "CalcPol Lorem Ipsum\n");
    python::def("CalcMR", calcMR,
        "CalcMR Lorem Ipsum\n");
    python::def("CalcFlexibility", calcFlexibility,
        "CalcFlexibility Lorem Ipsum\n");
    python::def("CalcODT", calcODT,
        "CalcODT Lorem Ipsum\n");
    python::def("CalcSchultz", calcSchultz,
        "CalcSchultz Lorem Ipsum\n");
    python::def("CalcRNCGRPCG", calcRNCG_RPCG,
        "CalcRNCG_RPCG Lorem Ipsum\n");
    python::def("CalcAZV", calcAZV,
        "CalcAZV Lorem Ipsum\n");
    python::def("CalcASV", calcASV,
        "CalcASV Lorem Ipsum\n");
    python::def("CalcDSV", calcDSV,
        "CalcDSV Lorem Ipsum\n");
    python::def("CalcAZS", calcAZS,
        "CalcAZS Lorem Ipsum\n");
    python::def("CalcASZ", calcASZ,
        "CalcASZ Lorem Ipsum\n");
    python::def("CalcDN2S", calcDN2S,
        "CalcDN2S Lorem Ipsum\n");
    python::def("CalcDN2I", calcDN2I,
        "CalcDN2I Lorem Ipsum\n");
    python::def("CalcASI", calcASI,
        "CalcASI Lorem Ipsum\n");
    python::def("CalcDSI", calcDSI,
        "CalcDSI Lorem Ipsum\n");
    python::def("CalcASN", calcASN,
        "CalcASN Lorem Ipsum\n");
    python::def("CalcDSN", calcDSN,
        "CalcDSN Lorem Ipsum\n");
    python::def("CalcDN2N", calcDN2N,
        "CalcDN2N Lorem Ipsum\n");
    python::def("CalcANS", calcANS,
        "CalcANS Lorem Ipsum\n");
    python::def("CalcANV", calcANV,
        "CalcANV Lorem Ipsum\n");
    python::def("CalcAZN", calcAZN,
        "CalcAZN Lorem Ipsum\n");
    python::def("CalcANZ", calcANZ,
        "CalcANZ Lorem Ipsum\n");
    python::def("CalcANI", calcANI,
        "CalcANI Lorem Ipsum\n");
    python::def("CalcDSZ", calcDSZ,
        "CalcDSZ Lorem Ipsum\n");
    python::def("CalcANN", calcANN,
        "CalcANN Lorem Ipsum\n");
    python::def("CalcDN2Z", calcDN2Z,
        "CalcDN2Z Lorem Ipsum\n");
    python::def("CalcANMat", calcANMat,
        "CalcANMat Lorem Ipsum\n");
    python::def("CalcAZMat", calcAZMat,
        "CalcAZMat Lorem Ipsum\n");
    python::def("CalcASMat", calcASMat,
        "CalcASMat Lorem Ipsum\n");
    python::def("CalcDSMat", calcDSMat,
        "CalcDSMat Lorem Ipsum\n");
    python::def("CalcDN2Mat", calcDN2Mat,
        "CalcDN2Mat Lorem Ipsum\n");
    python::def("CalcFrags", calcFrags,
        "CalcFrags Lorem Ipsum\n");
    python::def("CalcAddFeatures", calcAddFeatures,
        "CalcAddFeatures Lorem Ipsum\n");
    python::def("CalcInformationContent", calcInformationContent,
        "CalcInformationContent Lorem Ipsum\n");

};