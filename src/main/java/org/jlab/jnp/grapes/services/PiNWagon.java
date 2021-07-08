package org.jlab.jnp.grapes.services;

import org.jlab.jnp.physics.LorentzVector;

import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.data.Schema;
import org.jlab.jnp.hipo4.data.SchemaFactory;
import org.jlab.jnp.pdg.PDGDatabase;

import org.jlab.groot.data.H1F;

import java.util.Arrays;
import java.util.Map;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.IntStream;
import java.util.OptionalInt;

import org.jlab.jnp.utils.json.Json;
import org.jlab.jnp.utils.json.JsonObject;



/**
 * 
 * PiN Skimming
 *
 * @author Valerii
 */

public class PiNWagon extends HistoWagon {

	static final double PionMass   = 0.13957f;
	static final double ProtonMass = 0.938272;
   private double beamEnergy = 0;


   H1F hmm = new H1F("hmm","MM epi+X",200,0,2);

	public PiNWagon() {
		super("PiNWagon","valerii","0.0");
      register(hmm);
	}



   private static LorentzVector getLorentzVectorWithPID(Bank bnk, int ii, int pid) {
     float px = bnk.getFloat("px", ii);
     float py = bnk.getFloat("py", ii);
     float pz = bnk.getFloat("pz", ii);
     LorentzVector lvec = new LorentzVector();
     double mass = PDGDatabase.getParticleById(pid).mass();
     lvec.setPxPyPzM(px,py,pz,mass);
     return lvec;
   }



	@Override
	public boolean processDataEvent(Event event, SchemaFactory factory) {
		LorentzVector beam = new LorentzVector(0,0,beamEnergy,beamEnergy);
		LorentzVector targ = new LorentzVector(0,0,0,ProtonMass);

		Bank recPart = new Bank(factory.getSchema("REC::Particle"));
		Bank cnfPart = new Bank(factory.getSchema("RUN::config"));
		Bank calPart = new Bank(factory.getSchema("REC::Calorimeter"));
        
        Bank evPart = new Bank(factory.getSchema("REC::Event"));
        Bank sciPart = new Bank(factory.getSchema("REC::Scintillator"));
        Bank trajPart = new Bank(factory.getSchema("REC::Traj"));
        Bank tracPart = new Bank(factory.getSchema("REC::Track")); 
        
		event.read(recPart);
		event.read(cnfPart);
		event.read(calPart);
        
        event.read(evPart);
        event.read(sciPart);
        event.read(trajPart);
        event.read(tracPart);

		if( cnfPart!=null && recPart!=null && recPart.getRows()>1 ){
            int run = cnfPart.getInt("run",0);
            int evid = cnfPart.getInt("event", 0);

            List<Integer> ieles = new ArrayList<>();
            List<Integer> ipips = new ArrayList<>();
            
            bool only_epip = true;
             for(int ipart=0; ipart<recPart.getRows(); ipart++) {
              short stat = recPart.getShort("status", ipart);
              int pid = recPart.getInt("pid", ipart);
              int charge =  recPart.getInt("charge", ipart);
              if (charge != 0 && pid != 11 && pid != 211) only_epip = false;
             }
            
            if (!only_epip) continue;

            for(int ipart=0; ipart<recPart.getRows(); ipart++) {
              short stat = recPart.getShort("status", ipart);
              int pid = recPart.getInt("pid", ipart);
              if(stat<0 && pid==11) {
                ieles.add(ipart);
              } else if(pid==211 && stat>2000 && stat<4000) {
                  
                  
                double pipmom2 = Math.pow(recPart.getFloat("px",ipart), 2) + Math.pow(recPart.getFloat("py",ipart), 2) + Math.pow(recPart.getFloat("pz",ipart), 2);
                double pipmom = Math.sqrt(pipmom2);
                if (pipmom > 0.1) {
                  ipips.add(ipart);
                }
              }
            }

            BankBuilder builder = new BankBuilder(factory.getSchema("EXCLUSIVE::EPipN"));
            for(Integer iele: ieles) {
              for(Integer ipip: ipips) {
                OptionalInt esec = IntStream.range(0,calPart.getRows()).filter(ii -> calPart.getShort("pindex",ii)==iele).map(ii -> calPart.getByte("sector",ii)).findFirst();
                LorentzVector ele = getLorentzVectorWithPID(recPart,iele,11);
                LorentzVector pip = getLorentzVectorWithPID(recPart,ipip,211);
                LorentzVector epipX = beam.add(targ).sub(ele).sub(pip);
                if(epipX.mass()<4.) {
                  hmm.fill(epipX.mass());
                
                float eleDC_x[] = new float[3], eleDC_y[]  = new float[3], eleDC_z[]  = new float[3];
                float pipDC_x[] = new float[3], pipDC_y[] = new float[3], pipDC_z[] = new float[3];
                    
                bool elelayer_1 = false, elelayer_2 = false, elelayer_3 = false;
                bool piplayer_1 = false, piplayer_2 = false, piplayer_3 = false;
                    
                 for(int iTrajPart=0; iTrajPart<trajPart.getRows(); iTrajPart++) {
                     
                     int pind = trajPart.getInt("pindex", iTrajPart);
                     int detec = trajPart.getInt("detector", iTrajPart);
                     int layer = trajPart.getInt("layer", iTrajPart);
                     
                     if (detec == 6){
                         if (layer == 6){
                             if (pind == iele){
                                 elelayer_1 = true;
                                 eleDC_x[0] =  trajPart.getFloat("x", iTrajPart);
                                 eleDC_y[0] =  trajPart.getFloat("y", iTrajPart);
                                 eleDC_z[0] =  trajPart.getFloat("z", iTrajPart);
                             }
                             if (pind == ipip){
                                 piplayer_1 = true;
                                 pipDC_x[0] =  trajPart.getFloat("x", iTrajPart);
                                 pipDC_y[0] =  trajPart.getFloat("y", iTrajPart);
                                 pipDC_z[0] =  trajPart.getFloat("z", iTrajPart);
                             }
                         }
                         if (layer == 18){
                             if (pind == iele){
                                 elelayer_2 = true;
                                 eleDC_x[1] =  trajPart.getFloat("x", iTrajPart);
                                 eleDC_y[1] =  trajPart.getFloat("y", iTrajPart);
                                 eleDC_z[1] =  trajPart.getFloat("z", iTrajPart);
                             }
                             if (pind == ipip){
                                 piplayer_2 = true;
                                 pipDC_x[1] =  trajPart.getFloat("x", iTrajPart);
                                 pipDC_y[1] =  trajPart.getFloat("y", iTrajPart);
                                 pipDC_z[1] =  trajPart.getFloat("z", iTrajPart);
                             }
                         }
                         if (layer == 36){
                             if (pind == iele){
                                 elelayer_3 = true;
                                 eleDC_x[2] =  trajPart.getFloat("x", iTrajPart);
                                 eleDC_y[2] =  trajPart.getFloat("y", iTrajPart);
                                 eleDC_z[2] =  trajPart.getFloat("z", iTrajPart);
                             }
                             if (pind == ipip){
                                 piplayer_3 = true;
                                 pipDC_x[2] =  trajPart.getFloat("x", iTrajPart);
                                 pipDC_y[2] =  trajPart.getFloat("y", iTrajPart);
                                 pipDC_z[2] =  trajPart.getFloat("z", iTrajPart);
                             }
                         }
                         
                     }
                     
                 }
                    
                if (!(elelayer_1 * elelayer_2 * elelayer_3 * piplayer_1 * piplayer_2 * piplayer_3)) continue;

                  builder.addRow(new Object[][]{
                         {"ex", ele.px()},
                         {"ey", ele.py()},
                         {"ez", ele.pz()},
                         {"evz", recPart.getFloat("vz",iele)},
                         {"echi2pid", recPart.getFloat("chi2pid",iele)},
                         {"echelicity", evPart.getFloat("helicity",iele)},
                      
                         {"eDC_x_1", eleDC_x[0]},
                         {"eDC_x_2", eleDC_x[1]},
                         {"eDC_x_3", eleDC_x[2]},
                         {"eDC_y_1", eleDC_y[0]},
                         {"eDC_y_2", eleDC_y[1]},
                         {"eDC_y_3", eleDC_y[2]},
                         {"eDC_z_1", eleDC_z[0]},
                         {"eDC_z_2", eleDC_z[1]},
                         {"eDC_z_3", eleDC_z[2]},
                      
                         {"pipx", pip.px()},
                         {"pipy", pip.py()},
                         {"pipz", pip.pz()},
                         {"pipvz", recPart.getFloat("vz",ipip)},
                         {"pipchi2pid", recPart.getFloat("chi2pid",ipip)},
                      
                         {"pipDC_x_1", pipDC_x[0]},
                         {"pipDC_x_2", pipDC_x[1]},
                         {"pipDC_x_3", pipDC_x[2]},
                         {"pipDC_y_1", pipDC_y[0]},
                         {"pipDC_y_2", pipDC_y[1]},
                         {"pipDC_y_3", pipDC_y[2]},
                         {"pipDC_z_1", pipDC_z[0]},
                         {"pipDC_z_2", pipDC_z[1]},
                         {"pipDC_z_3", pipDC_z[2]},
                      
                         {"run", run},
                         {"evid", evid}
                  });
                }
              }
            }

            if(builder.validBank()) {
              event.write(builder.build());
              return true;
            }
          }

          return false;
     }
 

   
    @Override
    public final boolean init(String jsonString) {
        JsonObject jsonObj = Json.parse(jsonString).asObject();
        beamEnergy = jsonObj.getDouble("beamEnergy",-1.0);
        if (beamEnergy>0) {
          String beamSetup = "EB="+beamEnergy;
          System.out.println(engineName +" READY with "+beamSetup);
          return true;
        }
        System.out.println("Error initializing "+engineName+" due to beamEnergy YAML parameters.");
        return false;
    }

}
