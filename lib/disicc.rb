def import_my_data
  a1 = Alignment.get(2565)
  a2 = Alignment.get(2212)
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>130),
                              Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>367),
                              "temp_data/InterAlign/BEIV/BEIV.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>38),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>261),
    "temp_data/InterAlign/CDV/CDV.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>42),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>359),
    "temp_data/InterAlign/DMV/DMV.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>22),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>389),
    "temp_data/InterAlign/GPV/GPV.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>30),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>350),
    "temp_data/InterAlign/HEND/HEND.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>109),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>369),
    "temp_data/InterAlign/JV/JV.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>113),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>381),
    "temp_data/InterAlign/MENV/MENV.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>46),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>357),
    "temp_data/InterAlign/MeV/MeV.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>117),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>365),
    "temp_data/InterAlign/MOSV/MOSV.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>34),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>352),
    "temp_data/InterAlign/NIPH/NIPH.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>260),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>355),
    "temp_data/InterAlign/PDPR/PDPR.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>50),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>362),
    "temp_data/InterAlign/RPV/RPV.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>101),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>384),
    "temp_data/InterAlign/MuV/MuV.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>89),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>385),
    "temp_data/InterAlign/SPIV5/SPIV5.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>122),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>383),
    "temp_data/InterAlign/TIOV/TIOV.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>126),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>363),
    "temp_data/InterAlign/TUPV/TUPV.out")
end

def import_my_data1LP
  a1 = Alignment.get(2565)
  a2 = Alignment.get(2212)
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>130),
                              Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>367),
                              "temp_data/InterAlign/BEIV/BEIV.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>38),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>261),
    "temp_data/InterAlign/CDV/CDV.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>42),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>359),
    "temp_data/InterAlign/DMV/DMV.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>22),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>389),
    "temp_data/InterAlign/GPV/GPV.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>30),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>350),
    "temp_data/InterAlign/HEND/HEND.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>109),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>369),
    "temp_data/InterAlign/JV/JV.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>113),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>381),
    "temp_data/InterAlign/MENV/MENV.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>46),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>357),
    "temp_data/InterAlign/MeV/MeV.out")
end

def import_my_data2LP
  a1 = Alignment.get(2565)
  a2 = Alignment.get(2212)
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>117),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>365),
    "temp_data/InterAlign/MOSV/MOSV.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>34),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>352),
    "temp_data/InterAlign/NIPH/NIPH.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>260),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>355),
    "temp_data/InterAlign/PDPR/PDPR.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>50),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>362),
    "temp_data/InterAlign/PDV/PDV.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>54),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>358),
    "temp_data/InterAlign/RPV/RPV.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>101),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>384),
    "temp_data/InterAlign/MuV/MuV.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>89),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>385),
    "temp_data/InterAlign/SPIV5/SPIV5.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>122),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>383),
    "temp_data/InterAlign/TIOV/TIOV.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>126),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>363),
    "temp_data/InterAlign/TUPV/TUPV.out")
end

def run_inter
  
  puts "Starting Caps BEIV"
  system "cd lib/comp_apps/caps-perl/; ./caps.pl ../../../temp_data/InterAlign/BEIV/LNcaps.ctl"
  puts "Done Caps BEIV"
  puts "Starting Caps CDV"
  system "cd lib/comp_apps/caps-perl/; ./caps.pl ../../../temp_data/InterAlign/CDV/LNcaps.ctl"
  puts "Done Caps CDV"
  puts "Starting Caps DMV"
  system "cd lib/comp_apps/caps-perl/; ./caps.pl ../../../temp_data/InterAlign/DMV/LNcaps.ctl"
  puts "Done Caps DMV"
  puts "Starting Caps GPV"
  system "cd lib/comp_apps/caps-perl/; ./caps.pl ../../../temp_data/InterAlign/GPV/LNcaps.ctl"
  puts "Done Caps GPV"
  puts "Starting Caps HEND"
  system "cd lib/comp_apps/caps-perl/; ./caps.pl ../../../temp_data/InterAlign/HEND/LNcaps.ctl"
  puts "Done Caps HEND"
  puts "Starting Caps JV"
  system "cd lib/comp_apps/caps-perl/; ./caps.pl ../../../temp_data/InterAlign/JV/LNcaps.ctl"
  puts "Done Caps JV"
  puts "Starting Caps MENV"
  system "cd lib/comp_apps/caps-perl/; ./caps.pl ../../../temp_data/InterAlign/MENV/LNcaps.ctl"
  puts "Done Caps MENV"
  puts "Starting Caps MeV"
  system "cd lib/comp_apps/caps-perl/; ./caps.pl ../../../temp_data/InterAlign/MeV/LNcaps.ctl"
  puts "Done Caps MeV"
  puts "Starting Caps MOSV"
  system "cd lib/comp_apps/caps-perl/; ./caps.pl ../../../temp_data/InterAlign/MOSV/LNcaps.ctl"
  puts "Done Caps MOSV"
  puts "Starting Caps MuV"
  system "cd lib/comp_apps/caps-perl/; ./caps.pl ../../../temp_data/InterAlign/MuV/LNcaps.ctl"
  puts "Done Caps MuV"
  puts "Starting Caps NIPH"
  system "cd lib/comp_apps/caps-perl/; ./caps.pl ../../../temp_data/InterAlign/NIPH/LNcaps.ctl"
  puts "Done Caps NIPH"
  puts "Starting Caps PDV"
  system "cd lib/comp_apps/caps-perl/; ./caps.pl ../../../temp_data/InterAlign/PDV/LNcaps.ctl"
  puts "Done Caps PDV"
  puts "Starting Caps RPV"
  system "cd lib/comp_apps/caps-perl/; ./caps.pl ../../../temp_data/InterAlign/RPV/LNcaps.ctl"
  puts "Done Caps RPV"
  puts "Starting Caps SPIV5"
  system "cd lib/comp_apps/caps-perl/; ./caps.pl ../../../temp_data/InterAlign/SPIV5/LNcaps.ctl"
  puts "Done Caps SPIV5"
  puts "Starting Caps TIOV"
  system "cd lib/comp_apps/caps-perl/; ./caps.pl ../../../temp_data/InterAlign/TIOV/LNcaps.ctl"
  puts "Done Caps TIOV"
  puts "Starting Caps TUPV"
  system "cd lib/comp_apps/caps-perl/; ./caps.pl ../../../temp_data/InterAlign/TUPV/LNcaps.ctl"
  puts "Done Caps TUPV"
end

def run_inter1
  
  puts "Starting Caps BEIV"
  system "cd lib/comp_apps/caps-perl/; ./caps.pl ../../../temp_data/InterAlign/BEIV/LNcaps.ctl"
  puts "Done Caps BEIV"
  puts "Starting Caps CDV"
  system "cd lib/comp_apps/caps-perl/; ./caps.pl ../../../temp_data/InterAlign/CDV/LNcaps.ctl"
  puts "Done Caps CDV"
  puts "Starting Caps DMV"
  system "cd lib/comp_apps/caps-perl/; ./caps.pl ../../../temp_data/InterAlign/DMV/LNcaps.ctl"
  puts "Done Caps DMV"
  puts "Starting Caps GPV"
  system "cd lib/comp_apps/caps-perl/; ./caps.pl ../../../temp_data/InterAlign/GPV/LNcaps.ctl"
  puts "Done Caps GPV"
  puts "Starting Caps HEND"
  system "cd lib/comp_apps/caps-perl/; ./caps.pl ../../../temp_data/InterAlign/HEND/LNcaps.ctl"
  puts "Done Caps HEND"
  puts "Starting Caps JV"
  system "cd lib/comp_apps/caps-perl/; ./caps.pl ../../../temp_data/InterAlign/JV/LNcaps.ctl"
  puts "Done Caps JV"
  puts "Starting Caps MENV"
  system "cd lib/comp_apps/caps-perl/; ./caps.pl ../../../temp_data/InterAlign/MENV/LNcaps.ctl"
  puts "Done Caps MENV"
end

def run_inter2
  

  puts "Starting Caps MeV"
  system "cd lib/comp_apps/caps-perl/; ./caps.pl ../../../temp_data/InterAlign/MeV/LNcaps.ctl"
  puts "Done Caps MeV"
  puts "Starting Caps MOSV"
  system "cd lib/comp_apps/caps-perl/; ./caps.pl ../../../temp_data/InterAlign/MOSV/LNcaps.ctl"
  puts "Done Caps MOSV"
  puts "Starting Caps MuV"
  system "cd lib/comp_apps/caps-perl/; ./caps.pl ../../../temp_data/InterAlign/MuV/LNcaps.ctl"
  puts "Done Caps MuV"
  puts "Starting Caps NIPH"
  system "cd lib/comp_apps/caps-perl/; ./caps.pl ../../../temp_data/InterAlign/NIPH/LNcaps.ctl"
  puts "Done Caps NIPH"
  puts "Starting Caps PDV"
  system "cd lib/comp_apps/caps-perl/; ./caps.pl ../../../temp_data/InterAlign/PDV/LNcaps.ctl"
  puts "Done Caps PDV"
  puts "Starting Caps RPV"
  system "cd lib/comp_apps/caps-perl/; ./caps.pl ../../../temp_data/InterAlign/RPV/LNcaps.ctl"
  puts "Done Caps RPV"
  puts "Starting Caps SPIV5"
  system "cd lib/comp_apps/caps-perl/; ./caps.pl ../../../temp_data/InterAlign/SPIV5/LNcaps.ctl"
  puts "Done Caps SPIV5"
  puts "Starting Caps TIOV"
  system "cd lib/comp_apps/caps-perl/; ./caps.pl ../../../temp_data/InterAlign/TIOV/LNcaps.ctl"
  puts "Done Caps TIOV"
  puts "Starting Caps TUPV"
  system "cd lib/comp_apps/caps-perl/; ./caps.pl ../../../temp_data/InterAlign/TUPV/LNcaps.ctl"
  puts "Done Caps TUPV"
end
def import_my_data1LN
  a1 = Alignment.get(2565)
  a2 = Alignment.get(2804)
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>130),
                              Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>131),
                              "temp_data/InterAlign/BEIV/LN_BEIV.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>38),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>39),
    "temp_data/InterAlign/CDV/LN_CDV.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>42),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>43),
    "temp_data/InterAlign/DMV/LN_DMV.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>22),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>23),
    "temp_data/InterAlign/GPV/LN_GPV.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>30),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>31),
    "temp_data/InterAlign/HEND/LN_HEND.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>109),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>110),
    "temp_data/InterAlign/JV/LN_JV.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>113),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>114),
    "temp_data/InterAlign/MENV/LN_MENV.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>46),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>47),
    "temp_data/InterAlign/MeV/LN_MeV.out")
end

def import_my_data2LN
  a1 = Alignment.get(2565)
  a2 = Alignment.get(2804)
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>117),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>118),
    "temp_data/InterAlign/MOSV/LN_MOSV.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>34),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>35),
    "temp_data/InterAlign/NIPH/LN_NIPH.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>260),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>251),
    "temp_data/InterAlign/PDPR/LN_PDPR.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>50),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>51),
    "temp_data/InterAlign/PDV/LN_PDV.out")
    Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>54),
      Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>55),
      "temp_data/InterAlign/RPV/RPV.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>101),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>102),
    "temp_data/InterAlign/MuV/LN_MuV.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>89),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>93),
    "temp_data/InterAlign/SPIV5/LN_SPIV5.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>122),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>123),
    "temp_data/InterAlign/TIOV/LN_TIOV.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>126),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>127),
    "temp_data/InterAlign/TUPV/LN_TUPV.out")
end

def run_inter1PN
  
  puts "Starting Caps BEIV"
  system "cd lib/comp_apps/caps-perl/; ./caps.pl ../../../temp_data/InterAlign/BEIV/PNcaps.ctl"
  puts "Done Caps BEIV"
  puts "Starting Caps CDV"
  system "cd lib/comp_apps/caps-perl/; ./caps.pl ../../../temp_data/InterAlign/CDV/PNcaps.ctl"
  puts "Done Caps CDV"
  puts "Starting Caps DMV"
  system "cd lib/comp_apps/caps-perl/; ./caps.pl ../../../temp_data/InterAlign/DMV/PNcaps.ctl"
  puts "Done Caps DMV"
  puts "Starting Caps GPV"
  system "cd lib/comp_apps/caps-perl/; ./caps.pl ../../../temp_data/InterAlign/GPV/PNcaps.ctl"
  puts "Done Caps GPV"
  puts "Starting Caps HEND"
  system "cd lib/comp_apps/caps-perl/; ./caps.pl ../../../temp_data/InterAlign/HEND/PNcaps.ctl"
  puts "Done Caps HEND"
  puts "Starting Caps JV"
  system "cd lib/comp_apps/caps-perl/; ./caps.pl ../../../temp_data/InterAlign/JV/PNcaps.ctl"
  puts "Done Caps JV"
  puts "Starting Caps MENV"
  system "cd lib/comp_apps/caps-perl/; ./caps.pl ../../../temp_data/InterAlign/MENV/PNcaps.ctl"
  puts "Done Caps MENV"
end

def run_inter2PN
  

  puts "Starting Caps MeV"
  system "cd lib/comp_apps/caps-perl/; ./caps.pl ../../../temp_data/InterAlign/MeV/PNcaps.ctl"
  puts "Done Caps MeV"
  puts "Starting Caps MOSV"
  system "cd lib/comp_apps/caps-perl/; ./caps.pl ../../../temp_data/InterAlign/MOSV/PNcaps.ctl"
  puts "Done Caps MOSV"
  puts "Starting Caps MuV"
  system "cd lib/comp_apps/caps-perl/; ./caps.pl ../../../temp_data/InterAlign/MuV/PNcaps.ctl"
  puts "Done Caps MuV"
  puts "Starting Caps NIPH"
  system "cd lib/comp_apps/caps-perl/; ./caps.pl ../../../temp_data/InterAlign/NIPH/PNcaps.ctl"
  puts "Done Caps NIPH"
  puts "Starting Caps PDV"
  system "cd lib/comp_apps/caps-perl/; ./caps.pl ../../../temp_data/InterAlign/PDV/PNcaps.ctl"
  puts "Done Caps PDV"
  puts "Starting Caps RPV"
  system "cd lib/comp_apps/caps-perl/; ./caps.pl ../../../temp_data/InterAlign/RPV/PNcaps.ctl"
  puts "Done Caps RPV"
  puts "Starting Caps SPIV5"
  system "cd lib/comp_apps/caps-perl/; ./caps.pl ../../../temp_data/InterAlign/SPIV5/PNcaps.ctl"
  puts "Done Caps SPIV5"
  puts "Starting Caps TIOV"
  system "cd lib/comp_apps/caps-perl/; ./caps.pl ../../../temp_data/InterAlign/TIOV/PNcaps.ctl"
  puts "Done Caps TIOV"
  puts "Starting Caps TUPV"
  system "cd lib/comp_apps/caps-perl/; ./caps.pl ../../../temp_data/InterAlign/TUPV/PNcaps.ctl"
  puts "Done Caps TUPV"
end


def import_my_data1PN
  a1 = Alignment.get(2212)
  a2 = Alignment.get(2803)
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>367),
                              Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>131),
                              "temp_data/InterAlign/BEIV/PPNN_BEIV.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>261),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>39),
    "temp_data/InterAlign/CDV/PN_CDV.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>359),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>43),
    "temp_data/InterAlign/DMV/PN_DMV.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>389),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>23),
    "temp_data/InterAlign/GPV/PN_GPV.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>350),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>31),
    "temp_data/InterAlign/HEND/PN_HEND.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>369),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>110),
    "temp_data/InterAlign/JV/PN_JV.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>381),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>114),
    "temp_data/InterAlign/MENV/PN_MENV.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>357),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>47),
    "temp_data/InterAlign/MeV/PN_MeV.out")
end

def import_my_data2PN
  a1 = Alignment.get(2212)
  a2 = Alignment.get(2804)
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>365),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>118),
    "temp_data/InterAlign/MOSV/PN_MOSV.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>352),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>35),
    "temp_data/InterAlign/NIPH/PN_NIPH.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>355),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>251),
    "temp_data/InterAlign/PDPR/PN_PDPR.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>362),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>51),
    "temp_data/InterAlign/PDV/PN_PDV.out")
    Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>358),
      Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>55),
      "temp_data/InterAlign/RPV/RPV.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>384),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>102),
    "temp_data/InterAlign/MuV/PN_MuV.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>385),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>93),
    "temp_data/InterAlign/SPIV5/PN_SPIV5.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>383),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>123),
    "temp_data/InterAlign/TIOV/PN_TIOV.out")
  Alignment.import_inter_caps(Alignment.first(:alignment_name=>a1.alignment_name, :seq_id=>363),
    Alignment.first(:alignment_name=>a2.alignment_name, :seq_id=>127),
    "temp_data/InterAlign/TUPV/PN_TUPV.out")
end