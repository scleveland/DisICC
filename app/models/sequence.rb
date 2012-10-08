

class Sequence 
  include DataMapper::Resource
  
  property :seq_id, Serial
  property :seq_name, String, :required => true
  property :sequence, Text, :required => false
  property :seq_type, String, :required => true
  property :seq_accession, String, :required => true
  property :abrev_name, String, :required => false
  property :disorder_percent, Integer, :required => false
  property :alternate_name, String, :required => false
  
  
  #has n, :a_asequences
  #has n, :disorder
  
  def generate_aasequences
    (0..self.sequence.length-1).each do |i|
      AAsequence.create(:seq_id=> self.seq_id,
                       :amino_acid=>self.sequence[i],
                       :original_position=>i)
    end
  end
  
  
  def run_and_store_disorder
    begin
      self.generate_iupred_disorder_short
    rescue Exception => e  
      puts self.abrev_name + " IUPRED short failed!"
      puts e.message
    end
    begin
      self.generate_iupred_disorder
    rescue Exception => e  
      puts self.abrev_name + " IUPRED long failed!"
      puts e.message
    end
    begin
      self.generate_ronn
    rescue Exception => e  
      puts self.abrev_name + " RONN failed!"
      puts e.message
    end
    begin
      self.generate_pondr_fit
    rescue Exception => e  
      puts self.abrev_name + " PONDRFit failed!"
      puts e.message
    end
    begin
      self.generate_disembl
    rescue Exception => e  
      puts self.abrev_name + " DisEMBL failed!"
      puts e.message
    end
  end
  
  def generate_fasta_file
    filepath = "temp_data/"+self.abrev_name+"_"+self.seq_type+".fasta"
    f = File.new(filepath, "w+")
    f.write(">"+self.abrev_name+"|"+self.seq_name+"|"+self.seq_type+"|"+self.seq_accession+"\n")
    f.write(self.sequence)
    f.close
    return filepath
  end
  
  def generate_fasta_string
    fasta_string = ">"+self.abrev_name+"|"+self.seq_name+"|"+self.seq_type+"|"+self.seq_accession+"\n"
    fasta_string = fasta_string + self.sequence + "\n"
  end
  #### DISORDER ####
  
  def generate_iupred_disorder
    res = `./lib/disorder_apps/iupred/iupred_mac #{self.generate_fasta_file} long`
    filepath = "temp_data/"+self.abrev_name+"_"+self.seq_type+"_iupred.fasta"
    f = File.new(filepath, "w+")
    f.write(res)
    f.close
    self.store_iupred(filepath)
  end
  
  def generate_iupred_disorder_short
    res = `./lib/disorder_apps/iupred/iupred_mac #{self.generate_fasta_file} short`
    filepath = "temp_data/"+self.abrev_name+"_"+self.seq_type+"_iupred_short.fasta"
    f = File.new(filepath, "w+")
    f.write(res)
    f.close
    self.store_iupred_short(filepath)
  end
  
  def store_iupred(filepath)
    #create a new disorder object
    dis = Disorder.create(:seq_id => self.seq_id, :disorder_type=>"IUPred", :version=>1)
    file = File.new(filepath, 'r')
    counter = 1
    aa_count = 0
    while (line = file.gets)
      #puts "#{counter}: #{line}"
      counter = counter + 1
      if counter > 10
        line_array = line.split(' ')
        if aa = AAsequence.first(:seq_id => self.seq_id, :original_position=>aa_count, :amino_acid=>line_array[1])
        #puts "Amino Acid -"+line_array[1]+ " : " + aa.amino_acid + " | " + aa_count.to_s
        #if aa.amino_acid == line_array[1]  
          dv = DisorderValue.create(:disorder_id => dis.disorder_id, :aasequence_id => aa.AAsequence_id, :dvalue=>line_array[2].to_f) 
        end
        aa_count +=1
      end
    end
  end
  
  def store_iupred_short(filepath)
    #create a new disorder object
    dis = Disorder.create(:seq_id => self.seq_id, :disorder_type=>"IUPred Short", :version=>1)
    file = File.new(filepath, 'r')
    counter = 1
    aa_count = 0
    while (line = file.gets)
      #puts "#{counter}: #{line}"
      counter = counter + 1
      if counter > 10
        line_array = line.split(' ')
        if aa = AAsequence.first(:seq_id => self.seq_id, :original_position=>aa_count, :amino_acid=>line_array[1])
        #puts "Amino Acid -"+line_array[1]+ " : " + aa.amino_acid + " | " + aa_count.to_s
        #if aa.amino_acid == line_array[1]  
          dv = DisorderValue.create(:disorder_id => dis.disorder_id, :aasequence_id => aa.AAsequence_id, :dvalue=>line_array[2].to_f) 
        end
        aa_count +=1
      end
    end
  end
  
  def calculate_disorder_consensus(disorder_types)
    AAsequence.all(:seq_id => self.seq_id, :order =>[:original_position]).each do |aa|
      dis_sum = 0
      #disorder_types.each do |dis_type|
          dis_ids = Disorder.all(:disorder_type=>disorder_types, :seq_id=>self.seq_id).map{|k| k.disorder_id}
          dvs = DisorderValue.all(:aasequence_id => aa.AAsequence_id, :disorder_id=>dis_ids).map{|c| c.dvalue}
          dis_sum = dvs.sum
        if disorder_types.include?("DisEMBL Hotloops")
          dis_sum = dis_sum + 0.38 
          dis_num = dis_ids.length
        end
      #end
      aa.disorder_consensus = dis_sum/dis_num
      aa.save
    end
  end
  
  def self.calculate_disorder_consensus_for_types(ptype, disorder_types)
    Sequence.all(:seq_type=>ptype).each do |seq|
      seq.calculate_disorder_consensus(disorder_types)
    end
  end
  
  def generate_ronn
     url ="http://www.strubi.ox.ac.uk/RONN"
     cur = Curl::Easy.new(url)
     post_params= "sequence=#{'>'+self.abrev_name}\r\n#{self.sequence}&display_probs=y"
     cur.http_post(post_params)
     puts post_params
     puts cur.body_str.to_s
     s =cur.body_str.to_s.split('<pre>')
     a = s[2].split("</pre>")
     filepath= "temp_data/#{self.abrev_name}_RONN"
     f = File.new(filepath, "w+")
     f.write(a[0].to_s)
     f.close 
     puts filepath
     self.store_ronn(filepath)
   end

   def store_ronn(filepath)
     #create a new disorder object
     dis = Disorder.create(:seq_id => self.seq_id, :disorder_type=>"RONN", :version=>1)
     file = File.new(filepath, 'r')
     counter = 1
     aa_count = 0
     while (line = file.gets)
       #puts "#{counter}: #{line}"
       counter = counter + 1
       if counter > 2
         line_array = line.split
         if aa = AAsequence.first(:seq_id => self.seq_id, :original_position=>aa_count)
         #puts "Amino Acid -"+line_array[1]+ " : " + aa.amino_acid + " | " + aa_count.to_s
         #if aa.amino_acid == line_array[1]  
           DisorderValue.create(:disorder_id => dis.disorder_id, :aasequence_id => aa.AAsequence_id, :dvalue=>line_array[1].to_f) 
         end
         aa_count +=1
       end
     end
   end

   def self.generate_ronn_for_all(seq_type)
     self.all(:seq_type => seq_type).each do |seq|
       begin
         if Disorder.first(:seq_id=>seq.seq_id, :disorder_type=>"RONN")
           puts seq.abrev_name + " Already Exists: #{seq.seq_id}"
         else
           seq.generate_ronn
         end
       rescue
         puts seq.abrev_name + " Didn't Store: #{seq.seq_id}"
       end
     end
   end

   def generate_pondr_fit
     url ="http://www.disprot.org/action_predict.php"
     cur = Curl::Easy.new(url)
     post_params= "PONDRFIT=yes&native_sequence=#{'>'+self.abrev_name}\r\n#{self.sequence}&fontsize=small&plotwidth=7&xticincrement=100&plotheight=auto&filetype=eps&legend=full"
     cur.http_post(post_params)
     puts post_params
     puts cur.body_str.to_s
     s =cur.body_str.to_s.split('IUPRED SHORT DATA</a><br><a href=')
     f = s[1].split(">PONDR-FIT DATA")
     a = f[0].to_s.gsub('"','')
     file_url="http://www.disprot.org/" + a
     file_cur= Curl::Easy.new(file_url)
     file_cur.http_post()
     filepath= "temp_data/#{self.abrev_name}_PondrFit"
     f = File.new(filepath, "w+")
     f.write(file_cur.body_str)
     f.close 
     puts filepath
     self.store_pondr_fit(filepath)
   end

   def store_pondr_fit(filepath)
     #create a new disorder object
     dis = Disorder.create(:seq_id => self.seq_id, :disorder_type=>"PONDR Fit", :version=>1)
     file = File.new(filepath, 'r')
     counter = 1
     aa_count = 0
     while (line = file.gets)
       #puts "#{counter}: #{line}"
       counter = counter + 1
         line_array = line.split
         if aa = AAsequence.first(:seq_id => self.seq_id, :original_position=>aa_count)
         #puts "Amino Acid -"+line_array[1]+ " : " + aa.amino_acid + " | " + aa_count.to_s
         #if aa.amino_acid == line_array[1]  
           DisorderValue.create(:disorder_id => dis.disorder_id, :aasequence_id => aa.AAsequence_id, :dvalue=>line_array[2].to_f) 
         end
         aa_count +=1
     end
   end

   def self.generate_pondr_fit_for_all(seq_type)
     self.all(:seq_type => seq_type).each do |seq|
       begin
         if Disorder.first(:seq_id=>seq.seq_id, :disorder_type=>"PONDR Fit")
           puts seq.abrev_name + " Already Exists: #{seq.seq_id}"
         else
           seq.generate_pondr_fit
         end
       rescue
         puts seq.abrev_name + " Didn't Store: #{seq.seq_id}"
       end
     end
   end

   def generate_disembl
     url = "http://dis.embl.de/cgiDict.py"
     cur = Curl::Easy.new(url)
     post_params= "key=process&SP_entry=&sequence_string=>#{self.sequence.gsub("\r",'').gsub("\n",'')}&smooth_frame=8&peak_frame=8&join_frame=4&fold_coils=1.20&fold_rem465=1.20&fold_hotloops=1.40&plot_title=&tango_PH=7.40&tango_T=298.15&tango_I=0.02&tango_TFE=0.00&fold_tango=1.00"
     puts post_params
     cur.http_post(post_params)
     s = cur.body_str.to_s.split('<th>Download predictions</th>')
     f = s[1].split(">smoothed scores</a>")
     a = f[0].to_s.gsub("\n<td><a href=",'').gsub('"','')
     file_url = "http://dis.embl.de/" + a

     file_cur= Curl::Easy.new(file_url)
     file_cur.http_post()
     filepath= "temp_data/#{self.abrev_name}_Disembl"
     f = File.new(filepath, "w+")
     f.write(file_cur.body_str)
     f.close #this will be the smoothed scores file
     self.store_disembl(filepath)
   end

   def store_disembl(filepath)
     #create a new disorder object
     dis_hl = Disorder.create(:seq_id => self.seq_id, :disorder_type=>"DisEMBL Hotloops", :version=>1)
     dis_coil = Disorder.create(:seq_id => self.seq_id, :disorder_type=>"DisEMBL Coils", :version=>1)
     file = File.new(filepath, 'r')
     counter = 1
     aa_count = 0
     while (line = file.gets)
       #puts "#{counter}: #{line}"
       counter = counter + 1
       if counter > 2
         line_array = line.split('      ')
         if aa = AAsequence.first(:seq_id => self.seq_id, :original_position=>aa_count)
         #puts "Amino Acid -"+line_array[1]+ " : " + aa.amino_acid + " | " + aa_count.to_s
         #if aa.amino_acid == line_array[1]  
           DisorderValue.create(:disorder_id => dis_coil.disorder_id, :aasequence_id => aa.AAsequence_id, :dvalue=>line_array[0].to_f) 
           DisorderValue.create(:disorder_id => dis_hl.disorder_id, :aasequence_id => aa.AAsequence_id, :dvalue=>line_array[1].to_f) 
         end
         aa_count +=1
       end
     end
   end

   def self.generate_disembl_for_all(seq_type)
     self.all(:seq_type => seq_type).each do |seq|
       begin
         if Disorder.first(:seq_id=>seq.seq_id, :disorder_type=>["DisEMBL Hotloops","DisEMBL Coils"])
           puts seq.abrev_name + " Already Exists: #{seq.seq_id}"
         else
           seq.generate_disembl
         end
       rescue
         puts seq.abrev_name + " Didn't Store: #{seq.seq_id}"
       end
     end
   end
  
  # jalview_string = ""
  # low disorder  ffff00        
  # avg disorder  ffcc00        
  # medium disorder ff9900        
  # highly disordered ff6600        
  # extremely disordered  ff0000        
  # no disorder 0
  
  
  ### Jalview ANNOTATIONS
  
  def self.generate_all_sequence_javliew_annotation_iupred(ptype)
    require 'csv'
    CSV.open("temp_data/jalview_#{ptype}_#{Time.now}.gff", "wb", {:col_sep => "\t"}) do |csv|
    csv <<["possible disorder","ffffcc"]
    csv <<["low disorder","ffff00"]
    csv <<["avg disorder","ffcc00"]
    csv <<["medium disorder","ff9900"]
    csv <<["highly disordered","ff6600"]
    csv <<["extremely disordered","ff0000"]
    csv <<["no disorder","0"]
    Sequence.all(:seq_type => ptype).each do |seq|
      if dis = Disorder.first(:seq_id=>seq.seq_id, :disorder_type=>"IUPred")
      counter = 1
      DisorderValue.all(:disorder_id=>dis.disorder_id,:dvalue.gte => 0.4, :order=>[:disorder_value_id]).each do |dv|
        if dv.dvalue > 0.4 && dv.dvalue < 0.5
          feature_type = "possible disorder"
        elsif dv.dvalue > 0.5 && dv.dvalue < 0.6
          feature_type = "low disorder"
        elsif dv.dvalue >= 0.6 && dv.dvalue < 0.7
          feature_type = "avg disorder"
        elsif dv.dvalue >= 0.7 && dv.dvalue < 0.8
          feature_type = "medium disorder"
        elsif dv.dvalue >= 0.8 && dv.dvalue < 0.9
          feature_type = "highly disordered"
        elsif dv.dvalue > 0.9
          feature_type = "extremely disordered"
        else
          feature_type = "no disorder"
        end
        unless seq.alternate_name.nil?
          
          csv << ["None", "#{seq.alternate_name.split('/')[0]}", -1, "#{dv.aasequence.original_position}", "#{dv.aasequence.original_position}", "#{feature_type}"]
        end
        counter+=1
      end
      end
     end
    end
  end
  
  
  
  def self.generate_all_sequence_javliew_annotation(ptype,disorder_type)
    require 'csv'
    CSV.open("temp_data/jalview_#{ptype}_#{Time.now}.gff", "wb", {:col_sep => "\t"}) do |csv|
    csv <<["possible disorder","ffffcc"]
    csv <<["low disorder","ffff00"]
    csv <<["avg disorder","ffcc00"]
    csv <<["medium disorder","ff9900"]
    csv <<["highly disordered","ff6600"]
    csv <<["extremely disordered","ff0000"]
    csv <<["no disorder","0"]
    Sequence.all(:seq_type => ptype).each do |seq|
      if dis = Disorder.first(:seq_id=>seq.seq_id, :disorder_type=>disorder_type)
      counter = 1
      DisorderValue.all(:disorder_id=>dis.disorder_id,:dvalue.gte => 0.4, :order=>[:disorder_value_id]).each do |dv|
        if dv.dvalue > 0.4 && dv.dvalue < 0.5
          feature_type = "possible disorder"
        elsif dv.dvalue > 0.5 && dv.dvalue < 0.6
          feature_type = "low disorder"
        elsif dv.dvalue >= 0.6 && dv.dvalue < 0.7
          feature_type = "avg disorder"
        elsif dv.dvalue >= 0.7 && dv.dvalue < 0.8
          feature_type = "medium disorder"
        elsif dv.dvalue >= 0.8 && dv.dvalue < 0.9
          feature_type = "highly disordered"
        elsif dv.dvalue > 0.9
          feature_type = "extremely disordered"
        else
          feature_type = "no disorder"
        end
        unless seq.alternate_name.nil?
          
          csv << ["None", "#{seq.alternate_name.split('/')[0]}", -1, "#{dv.aasequence.original_position}", "#{dv.aasequence.original_position}", "#{feature_type}"]
        end
        counter+=1
      end
      end
     end
    end
  end
  
  def self.generate_all_sequence_disorder_consensus_javliew_annotation(ptype,abrev_name=true)
    require 'csv'
    CSV.open("temp_data/jalview_#{ptype}_#{Time.now}.gff", "wb", {:col_sep => "\t"}) do |csv|
      csv <<["possible disorder","ffffcc"]
      csv <<["low disorder","ffff00"]
      csv <<["avg disorder","ffcc00"]
      csv <<["medium disorder","ff9900"]
      csv <<["highly disordered","ff6600"]
      csv <<["extremely disordered","ff0000"]
      csv <<["no disorder","0"]
      Sequence.all(:seq_type => ptype).each do |seq|
        AAsequence.all(:seq_id => seq.seq_id, :order=>[:original_position]).each do |aa|
          if aa.disorder_consensus > 0.4 && aa.disorder_consensus < 0.5
            feature_type = "possible disorder"
          elsif aa.disorder_consensus > 0.5 && aa.disorder_consensus < 0.6
            feature_type = "low disorder"
          elsif aa.disorder_consensus >= 0.6 && aa.disorder_consensus < 0.7
            feature_type = "avg disorder"
          elsif aa.disorder_consensus >= 0.7 && aa.disorder_consensus < 0.8
            feature_type = "medium disorder"
          elsif aa.disorder_consensus >= 0.8 && aa.disorder_consensus < 0.9
            feature_type = "highly disordered"
          elsif aa.disorder_consensus > 0.9
            feature_type = "extremely disordered"
          else
            feature_type = "no disorder"
          end
          if abrev_name
            csv << ["None", "#{seq.abrev_name}", -1, "#{aa.original_position}", "#{aa.original_position}", "#{feature_type}"]
          else          
            csv << ["None", "#{seq.alternate_name.split('/')[0]}", -1, "#{aa.original_position}", "#{aa.original_position}", "#{feature_type}"]
          
          end
        end
      end
    end
  end
  
  def generate_jalview_annotation_iupred
    jalview_string= ""
    dis = Disorder.first(:seq_id=>self.seq_id)
    counter = 1
    DisorderValue.all(:disorder_id=>dis.disorder_id,:dvalue.gte => 0.5, :order=>[:disorder_value_id]).each do |dv|
      if dv.dvalue > 0.5 && dv.dvalue < 0.6
        feature_type = "low disorder"
      elsif dv.dvalue >= 0.6 && dv.dvalue < 0.7
        feature_type = "avg disorder"
      elsif dv.dvalue >= 0.7 && dv.dvalue < 0.8
        feature_type = "medium disorder"
      elsif dv.dvalue >= 0.8 && dv.dvalue < 0.9
        feature_type = "highly disordered"
      elsif dv.dvalue > 0.9
        feature_type = "extremely disordered"
      else
        feature_type = "no disorder"
      end
      jalview_string = jalview_string + "None #{self.abrev_name} -1 #{counter} #{counter} #{feature_type}"+"\n"
   end
   return jalview_string
  end
  
  def generate_jalview_annotation(disorder_type)
    jalview_string= ""
    dis = Disorder.first(:seq_id=>self.seq_id, :disorder_type => disorder_type)
    counter = 1
    DisorderValue.all(:disorder_id=>dis.disorder_id,:dvalue.gte => 0.5, :order=>[:disorder_value_id]).each do |dv|
      if dv.dvalue > 0.5 && dv.dvalue < 0.6
        feature_type = "low disorder"
      elsif dv.dvalue >= 0.6 && dv.dvalue < 0.7
        feature_type = "avg disorder"
      elsif dv.dvalue >= 0.7 && dv.dvalue < 0.8
        feature_type = "medium disorder"
      elsif dv.dvalue >= 0.8 && dv.dvalue < 0.9
        feature_type = "highly disordered"
      elsif dv.dvalue > 0.9
        feature_type = "extremely disordered"
      else
        feature_type = "no disorder"
      end
      jalview_string = jalview_string + "None #{self.abrev_name} -1 #{counter} #{counter} #{feature_type}"+"\n"
   end
   return jalview_string
  end
  
  def set_aasequence_defaults
      sql="UPDATE a_asequences SET disorder_consensus =0.0 WHERE seq_id=#{self.seq_id}"
      repository.adapter.execute(sql)
      sql="UPDATE a_asequences SET contact_consensus =0.0 WHERE seq_id=#{self.seq_id}"
      repository.adapter.execute(sql)
      sql="UPDATE a_asequences SET contact_positive_consensus =0.0 WHERE seq_id=#{self.seq_id}"
      repository.adapter.execute(sql)
  end
end

