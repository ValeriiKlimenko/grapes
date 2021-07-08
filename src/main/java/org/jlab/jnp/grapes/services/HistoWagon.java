package org.jlab.jnp.grapes.services;

import org.jlab.groot.data.TDirectory;
import org.jlab.groot.data.IDataSet;
import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Schema;
import java.lang.Number;
import java.util.List;
import java.util.ArrayList;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

/**
 *
 * 
 * @author kenjo
 */

public abstract class HistoWagon extends Wagon {
    protected Map<String, IDataSet> histos = new ConcurrentHashMap<>();

    public HistoWagon(String name, String author, String version){        
        super(name,author,version);
    }

    protected void register(IDataSet... dss) {
      for(IDataSet ds: dss) {
        histos.put(ds.getName(), ds);
      }
    }

    public static class BankBuilder {
      private final List<Object[][]> rows = new ArrayList<>();
      private final Schema schema;

      public BankBuilder(Schema schema) {
        this.schema = schema;
      }

      public BankBuilder addRow(Object[][] vars) {
        rows.add(vars);
        return this;
      }

      public Boolean validBank() {
        return rows.size()>0;
      }

      public Bank build() {
        Bank bank = new Bank(schema, rows.size());
        for(int irow=0;irow<rows.size();irow++) {
          for(int ivar=0; ivar<rows.get(irow).length; ivar++) {
            bank.putFloat((String) rows.get(irow)[ivar][0], irow, ((Number) rows.get(irow)[ivar][1]).floatValue());
          }
        }
        return bank;
      }
    }

    @Override
    public void destroy() {
        TDirectory out = new TDirectory();
        out.mkdir("/");
        out.cd("/");
        for(String name: histos.keySet()) {
          out.addDataSet(histos.get(name));
        }
        out.writeFile(engineName+".hipo");
    }
}
