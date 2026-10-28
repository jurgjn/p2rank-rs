import java.text.*;
import java.util.*;

public class FmtTest {
    public static void main(String[] a) throws Exception {
        DecimalFormatSymbols sym = new DecimalFormatSymbols(); sym.setDecimalSeparator('.');
        DecimalFormat f2 = new DecimalFormat("0.00", sym), f3 = new DecimalFormat("0.000", sym), f4 = new DecimalFormat("0.0000", sym);
        Random r = new Random(42);
        List<Double> vals = new ArrayList<>();
        double[] special = {0.0, -0.0, 0.125, 0.375, 1.0005, 2.5, -2.5, 0.005, 0.0005, 0.00005, 1e-7, -1e-7, 123456.789, 1e10, 0.1, 0.7, 9.995, 9.9995, -0.004, -0.0049, 0.0045, 1.0E-4};
        for (double d : special) vals.add(d);
        for (int i = 0; i < 200000; i++) {
            int kind = i % 6;
            double v;
            if (kind == 0) v = (r.nextDouble() - 0.5) * 200;
            else if (kind == 1) v = r.nextDouble();
            else if (kind == 2) v = Math.round(r.nextDouble() * 1e6) / 1e6 * (r.nextBoolean() ? 1 : -1) * 100;
            else if (kind == 3) v = (r.nextInt(200000) - 100000) / 2000.0;   // many exact ties in binary/decimal
            else if (kind == 4) v = (r.nextInt(2000000) - 1000000) / 16384.0;
            else v = Double.longBitsToDouble(r.nextLong() & 0x7fffffffffffffffL) % 1000;
            if (Double.isNaN(v) || Double.isInfinite(v)) continue;
            vals.add(v);
        }
        StringBuilder sb = new StringBuilder();
        for (double v : vals) {
            sb.setLength(0);
            sb.append(Long.toHexString(Double.doubleToRawLongBits(v))).append('|')
              .append(f2.format(v)).append('|').append(f3.format(v)).append('|').append(f4.format(v)).append('|')
              .append(String.format("%.4f", v)).append('|').append(String.format("%8.3f", v)).append('|')
              .append(String.format("%6.3f", v)).append('|').append(String.format("%.3f", v)).append('|')
              .append(Double.toString(v));
            System.out.println(sb);
        }
    }
}
