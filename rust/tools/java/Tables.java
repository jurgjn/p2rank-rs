import java.util.*;
import javax.vecmath.Point3d;
import org.openscience.cdk.geometry.surface.Tessellate;
import org.openscience.cdk.config.Elements;
import org.biojava.nbio.structure.Element;

public class Tables {
    static String h(double d) { return Long.toHexString(Double.doubleToRawLongBits(d)); }
    public static void main(String[] a) throws Exception {
        String what = a[0];
        if (what.equals("tess")) {
            int level = Integer.parseInt(a[1]);
            Tessellate t = new Tessellate("ico", level);
            t.doTessellate();
            Point3d[] pts = t.getTessAsPoint3ds();
            System.out.println("# level " + level + " n=" + pts.length + " triangles=" + t.getNumberOfTriangles());
            for (Point3d p : pts) System.out.println(h(p.x) + " " + h(p.y) + " " + h(p.z));
        } else if (what.equals("elements")) {
            // biojava element enum names; and for each, the CDK symbol validity + vdw radius
            java.lang.reflect.Method m = Class.forName("cz.cuni.cusbg.surface.FasterNumericalSurface").getDeclaredMethod("getVdwRadius", org.openscience.cdk.interfaces.IAtom.class);
            m.setAccessible(true);
            for (Element e : Element.values()) {
                String sym = e.name();
                boolean known = !Elements.ofString(sym).equals(Elements.Unknown);
                String r = "-";
                if (known) {
                    r = h((double) m.invoke(null, new org.openscience.cdk.silent.Atom(sym)));
                }
                System.out.println("E " + sym + " " + known + " " + r);
            }
            for (char c = ' '; c <= '~'; c++) {
                String sym = String.valueOf(c);
                boolean known = !Elements.ofString(sym).equals(Elements.Unknown);
                String r = "-";
                if (known) { try { r = h((double) m.invoke(null, new org.openscience.cdk.silent.Atom(sym))); } catch (Throwable t) { r = "THROWS"; } }
                System.out.println("C " + (int) c + " " + known + " " + r);
            }
        }
    }
}
