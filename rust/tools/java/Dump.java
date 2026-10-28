import java.io.*;
import java.util.*;
import cz.siret.prank.program.params.*;
import cz.siret.prank.domain.*;
import cz.siret.prank.domain.loaders.LoaderParams;
import cz.siret.prank.features.*;
import cz.siret.prank.geom.*;
import cz.siret.prank.program.ml.Model;
import cz.siret.prank.prediction.pockets.rescorers.InstancePredictor;
import org.biojava.nbio.structure.Atom;

/**
 * Dumps intermediate stages of P2Rank prediction with exact double bits.
 * usage: Dump <distroDir> <structureFile> <outPrefix>
 */
public class Dump {
    static String h(double d) { return Long.toHexString(Double.doubleToRawLongBits(d)); }

    public static void main(String[] args) throws Exception {
        String distro = args[0];
        String file = args[1];
        String out = args[2];

        Params p = Params.getInst();
        p.setInstallDir(distro);
        ConfigLoader.overrideConfig(p, new File(distro + "/config/default.groovy"));
        AminoAcidMapper.initialize(p.getAa_mapping());
        LoaderParams.ignoreLigandsSwitch = true;

        Protein prot = Protein.load(file, new LoaderParams());

        try (PrintWriter w = new PrintWriter(new FileWriter(out + ".atoms"))) {
            for (Atom a : prot.getProteinAtoms()) {
                w.println(a.getPDBserial() + " " + a.getName() + " " + a.getGroup().getPDBName() + " "
                        + a.getGroup().getResidueNumber().printFull() + " " + a.getElement() + " "
                        + h(a.getX()) + " " + h(a.getY()) + " " + h(a.getZ()) + " " + h(a.getTempFactor()));
            }
        }

        prot.calcuateSurfaceAndExposedAtoms();
        Atoms sas = prot.getAccessibleSurface().getPoints();
        try (PrintWriter w = new PrintWriter(new FileWriter(out + ".sas"))) {
            for (Atom a : sas) {
                w.println(h(a.getX()) + " " + h(a.getY()) + " " + h(a.getZ()));
            }
        }
        try (PrintWriter w = new PrintWriter(new FileWriter(out + ".exposed"))) {
            for (Atom a : prot.getExposedAtoms()) {
                w.println(a.getPDBserial());
            }
        }

        PrankFeatureExtractor factory = (PrankFeatureExtractor) FeatureExtractor.createFactory();
        PrankFeatureExtractor proto = (PrankFeatureExtractor) factory.createPrototypeForProtein(prot, null);
        try (PrintWriter w = new PrintWriter(new FileWriter(out + ".deep"))) {
            for (Atom a : proto.getDeepLayer()) {
                w.println(a.getPDBserial());
            }
        }
        FeatureExtractor ext = proto.createInstanceForWholeProtein(null);
        List<FeatureVector> vectors = new ArrayList<>();
        try (PrintWriter w = new PrintWriter(new FileWriter(out + ".vectors"))) {
            for (Atom pt : ext.getSampledPoints().getPoints()) {
                FeatureVector v = ext.calcFeatureVector(pt);
                vectors.add(v);
                StringBuilder sb = new StringBuilder();
                for (double d : v.getArray()) {
                    sb.append(h(d)).append(' ');
                }
                w.println(sb.toString().trim());
            }
        }

        if (args.length > 3) {
            Model model = Model.load(distro + "/models/default");
            InstancePredictor pred = InstancePredictor.create(model, proto);
            double[] scores = pred.predictBatch(vectors);
            try (PrintWriter w = new PrintWriter(new FileWriter(out + ".scores"))) {
                for (double d : scores) {
                    w.println(h(d));
                }
            }
        }
    }
}
