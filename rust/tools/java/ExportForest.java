import java.io.*;
import java.lang.reflect.Field;
import java.nio.*;
import java.nio.channels.FileChannel;
import java.nio.file.*;

public class ExportForest {
    static Object get(Object o, Class<?> c, String name) throws Exception {
        Field f = c.getDeclaredField(name); f.setAccessible(true); return f.get(o);
    }
    public static void main(String[] a) throws Exception {
        Object c = cz.siret.prank.utils.WekaUtils.loadClassifier(cz.siret.prank.utils.Futils.inputStream(a[0]));
        Class<?> fb = Class.forName("cz.siret.prank.fforest.api.FlatBinaryForest");
        Class<?> lfb = Class.forName("cz.siret.prank.fforest.api.LegacyFlatBinaryForest");
        if (!lfb.isInstance(c)) throw new RuntimeException("not LegacyFlatBinaryForest: " + c.getClass());
        int numTrees = (Integer) get(c, fb, "numTrees");
        int numAttributes = (Integer) get(c, fb, "numAttributes");
        int[] childLeft = (int[]) get(c, fb, "childLeft");
        int[] childRight = (int[]) get(c, fb, "childRight");
        int[] attributeIndex = (int[]) get(c, fb, "attributeIndex");
        double[] splitPoint = (double[]) get(c, fb, "splitPoint");
        double[][] classProbs = (double[][]) get(c, lfb, "classProbs");
        int n = childLeft.length;
        System.out.println("trees=" + numTrees + " attrs=" + numAttributes + " nodes=" + n + " leaves=" + classProbs.length
            + " cr=" + childRight.length + " ai=" + attributeIndex.length + " sp=" + splitPoint.length);
        for (int i = 0; i < classProbs.length; i++) {
            if (classProbs[i] != null && classProbs[i].length != 2) throw new RuntimeException("leaf " + i + " len " + classProbs[i].length);
        }
        long size = 8 + 4*4 + (long)n*4*3 + (long)n*8 + (long)classProbs.length*16;
        ByteBuffer b = ByteBuffer.allocate((int) size).order(ByteOrder.LITTLE_ENDIAN);
        b.put("P2RFLB01".getBytes());
        b.putInt(numTrees); b.putInt(numAttributes); b.putInt(n); b.putInt(classProbs.length);
        for (int v : childLeft) b.putInt(v);
        for (int v : childRight) b.putInt(v);
        for (int v : attributeIndex) b.putInt(v);
        for (double v : splitPoint) b.putDouble(v);
        for (double[] p : classProbs) {
            if (p == null) { b.putDouble(Double.NaN); b.putDouble(Double.NaN); }
            else { b.putDouble(p[0]); b.putDouble(p[1]); }
        }
        Files.write(Paths.get(a[1]), b.array());
    }
}
