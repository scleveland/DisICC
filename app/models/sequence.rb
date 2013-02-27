

class Sequence 
  include DataMapper::Resource
  
  property :seq_id, Serial#, :field=>'seq_id'
  property :seq_name, String, :required => true, :length => 256
  property :sequence, Text, :required => true
  property :seq_type, String, :required => false
  property :seq_accession, String, :required => false, :default=>"none"
  property :abrev_name, String, :required => false
  property :disorder_percent, Integer, :required => false
  property :alternate_name, String, :required => false
  property :owner, Integer, :required => true, :default => 1
  property :deleted_at, ParanoidDateTime
  
 # alias :seq_id :id
  alias :id :seq_id
  has n, :a_asequences, 'AAsequence', :child_key=>[:seq_id]
  has n, :users, :through =>Resource
  has n, :disorders, 'Disorder', :child_key=>[:seq_id]
  has n, :intra_residue_contacts, 'IntraResidueContact', :child_key=>[:seq_id]
  has n, :caps, 'Caps', :child_key=>[:seq_id]
  has n, :new_caps, 'NewCap', :child_key=>[:seq_id]
  has n, :xdets, :through => :a_asequences
  has n, :conseqs, :through => :a_asequences
  has n, :alignments, 'Alignment', :child_key=>[:seq_id]
  
  after :save do
    u = User.get(self.owner)
    unless  u.sequences.first(:seq_id => self.seq_id)
      u.sequences << self
      u.save
    end
  end
  
  def clean_up_abrev
    if (self.abrev_name.count "/") > 0
      puts "yes"
      self.abrev_name = self.abrev_name.split("/")[0]
      self.save
    end
  end
  
  # create_sequence_from_fasta
  # This reads a fasta file and stores the sequences in the database if it doesn't already 
  # exist.  This means that the sequence and definition must match exactly otherwise a new
  # instance will be stored.
  #
  # @example Sequence.create_sequences_from_fasta("temp/myFastaFile.fa")
  #
  # @param[String] fasta_filename, a path to avalid fasta file with amino acid sequences
  # 
  # @returns[Array] returns and array of DataMapper sequence objects
  # @author Sean Cleveland
  # @updated 2012/10/9
  def self.create_sequences_from_fasta(fasta_filename, owner=1)
    sequences = []
    file = File.new(fasta_filename, 'r')
    ff = Bio::FlatFile.new(Bio::FastaFormat, file)
    ff.each_entry do |f|
      puts f.definition
      sequences << Sequence.create_sequence_from_bioruby_fasta_entry(f, owner=1)
    end
    sequences
  end
  
  def self.create_sequence_from_bioruby_fasta_entry(fasta_entry, seq_type, owner=1)
    seq = self.new(:seq_name=>fasta_entry.definition, 
             :abrev_name=>fasta_entry.definition[0..3],
             :seq_type=>seq_type,
             :sequence => fasta_entry.naseq.gsub('-',''), 
             :owner=>owner,
             )
    seq.valid?
    puts seq.errors.inspect()
    seq.save
    seq.generate_aasequences()
    seq
  end
  
  def calculate_contact_consensus(thread_num=100)
    aa_array =[]
    self.a_asequences.each do |amino_acid|
      aa_array << amino_acid
    end
    thread_array=[]
    thread_num.times do |i|
        thread_array[i] = Thread.new{
          while aa_array.length > 0 do
            aaseq = aa_array.pop
            count = 0
            # if !IntraResidueContact.first(:seq_id => aaseq.seq_id, :first_residue=> aaseq.original_position).nil?
            #                 count +=1
            #               elsif !IntraResidueContact.first(:seq_id => aaseq.seq_id, :second_residue=> aaseq.original_position).nil?
            #                 count +=1
            #               end
            #               if !Conseq.first(:aasequence_id => aaseq.AAsequence_id).nil?
            #                 if Conseq.first(:aasequence_id => aaseq.AAsequence_id).color < 4
            #                   count +=1
            #                 end
            #               end
            if !Xdet.first(:aasequence_id => aaseq.AAsequence_id).nil?
              if Xdet.first(:aasequence_id => aaseq.AAsequence_id).correlation > 0.0 || Xdet.first(:aasequence_id => aaseq.AAsequence_id).correlation == -2
                count +=1
                puts "Xdet"
              end
            end
            if !NewCap.first(:seq_id=> aaseq.seq_id, :position_one => aaseq.original_position).nil?
              count +=1
            elsif !NewCap.first(:seq_id=> aaseq.seq_id, :position_two => aaseq.original_position).nil?
              count +=1
            end
            aaseq.contact_consensus = count /2#4
            puts count.to_s + " : " + aaseq.contact_consensus.to_s
            aaseq.save
          end  
        }
     end
     thread_array.map{|t| t.join}          
  end
  
  def calculate_intra_consensus(special=false, thread_num=100)
    aaseq_array = self.a_asequences.to_a 
    thread_array=[]
    thread_num.times do |i|
      thread_array[i] = Thread.new{
         while aaseq_array.length > 0 do
            aaseq = aaseq_array.pop
            if special
              aaseq.calculate_intra_consensus_value_special
            else
              aaseq.calculate_intra_consensus_value
            end
         end
      }
    end
    thread_array.map{|t| t.join}
  end
  
  def generate_aasequences
    self.a_asequences.destroy!
    sequence = self.sequence.gsub("\n","").gsub("\r","")
    (0..sequence.length-1).each do |i|
      AAsequence.create(:seq_id=> self.seq_id,
                       :amino_acid=>sequence[i],
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
  
  def generate_fasta_file_one_line
    filepath = "temp_data/"+self.abrev_name+"_"+self.seq_type+".fasta"
    f = File.new(filepath, "w+")
    f.write(">"+self.abrev_name + "\n")#"|"+self.seq_name+"|"+self.seq_type+"|"+self.seq_accession+"\n")
    #f.write(self.a_asequences(:order=>[:seq_id], :fields=>[:amino_acid]).map{|aa| aa.amino_acid}.join(''))
    f.write(self.sequence + "\n")
    f.close
    return filepath
  end
  
  def run_svmcon
    path = self.generate_fasta_file_one_line
    puts "Starting SVMCon #{self.abrev_name}"
    system "~/svmcon1.0/bin/predict_map.sh #{path} #{path}.map"
    puts "Finsihed SVMCon for #{self.abrev_name}"
  end
  
  def run_nncon(out_dir)
    path = self.generate_fasta_file_one_line
    puts "Starting NNCon #{self.abrev_name}"
    puts "~/nncon1.0/bin/predict_ss_sa_cm.sh #{path} #{out_dir}"
    system "~/nncon1.0/bin/predict_ss_sa_cm.sh #{path} #{out_dir}"
    puts "Finsihed NNCon for #{self.abrev_name}"
  end
  
  def generate_fasta_string
    fasta_string = ">"+self.abrev_name+"|"+self.seq_name+"|"+self.seq_type+"|"+self.seq_accession+"\n"
    fasta_string = fasta_string + self.sequence + "\n"
  end
  
  def generate_fasta_string_one_line
    fasta_string = ">"+self.abrev_name+"|"+self.seq_name+"|"+self.seq_type+"|"+self.seq_accession+"\n"
    fasta_string = fasta_string + self.a_asequences(:order=>[:seq_id], :fields=>[:amino_acid]).map{|aa| aa.amino_acid}.join('') + "\n"
  end
  
  def sequence_one_line
    fasta_string = self.a_asequences(:order=>[:seq_id], :fields=>[:amino_acid]).map{|aa| aa.amino_acid}.join('')
  end
  #### DISORDER ####
  
  def generate_iupred_disorder
    unless self.disorders.first(:disorder_type=>"IUPred")
      res = `./lib/disorder_apps/iupred/iupred_mac #{self.generate_fasta_file} long`
      filepath = "temp_data/"+self.abrev_name+"_"+self.seq_type+"_iupred.fasta"
      f = File.new(filepath, "w+")
      f.write(res)
      f.close
      self.store_iupred(filepath)
    else
      puts "IUPred Already Stored!***********************"
    end
  end
  
  def generate_iupred_disorder_short
    unless self.disorders.first(:disorder_type=>"IUPred Short")
      res = `./lib/disorder_apps/iupred/iupred_mac #{self.generate_fasta_file} short`
      filepath = "temp_data/"+self.abrev_name+"_"+self.seq_type+"_iupred_short.fasta"
      f = File.new(filepath, "w+")
      f.write(res)
      f.close
      self.store_iupred_short(filepath)
    else
      puts "IUPred Short Already Stored!***********************"
    end
  end
  
  def store_iupred(filepath)
    #create a new disorder object
    dis = Disorder.create(:seq_id => self.seq_id, :disorder_type=>"IUPred", :version=>1)
    self.disorders << dis
    self.save
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
    self.disorders << dis
    self.save
    file = File.new(filepath, 'r')
    counter = 1
    aa_count = 0
    while (line = file.gets)
      puts "#{counter}: #{line}"
      counter = counter + 1
      puts counter
      if counter > 10
        line_array = line.split(' ')
        if aa = AAsequence.first(:seq_id => self.seq_id, :original_position=>aa_count, :amino_acid=>line_array[1])
        puts "Amino Acid -"+line_array[1]+ " : " + aa.amino_acid + " | " + aa_count.to_s + "|" + line_array[2].to_f.to_s
        #if aa.amino_acid == line_array[1]  
          dv = DisorderValue.new(:disorder_id => dis.disorder_id, :aasequence_id => aa.AAsequence_id, :dvalue=>line_array[2].to_f) 
          dv.valid?
          puts dv.errors.inspect()
          dv.save
        end
        aa_count +=1
      end
    end
  end
  
  #  def calculate_disorder_consensus()#disorder_types)
  #   dis_ids = self.disorders.map{|k| k.disorder_id}
  #   puts dis_ids.join(',')
  #   AAsequence.all(:seq_id => self.seq_id, :order =>[:original_position]).each do |aa|
  #     # dis_sum = 0
  #     #       dvals = aa.disorder_values#DisorderValue.all(:aasequence_id => aa.id, :disorder_id=>dis_ids)
  #     #       dvs = dvals.map{|c| c.dvalue}
  #     #       puts "DVals: " + dvs.join(',').to_s
  #     #       dis_sum = dvs.inject{|sum,x| sum + x }
  #     #       puts "Sum: " + dis_sum.to_s
  #     #         # if disorder_types.include?("DisEMBL Hotloops")
  #     #         #   dis_sum = dis_sum + 0.38 
  #     #         #   dis_num = dis_ids.length
  #     #         # end
  #     #       #end
  #     #       if dvals.count > 0
  #     #         aa.disorder_consensus = dis_sum/dvals.count
  #     #       else
  #     #         aa.disorder_consensus = 0
  #     #       end
  #     aa.disorder_consensus = aa.disorder_values.avg(:dvalue)
  #     puts "Consensus: "+aa.disorder_consensus.to_s
  #     aa.save
  #   end
  # end
  
  def calculate_disorder_consensus_threaded(thread_num = 100)#disorder_types)
    aa_array =[]
    AAsequence.all(:seq_id => self.seq_id, :order =>[:original_position]).each do |aa|
      aa_array << aa
    end
    dis_ids = self.disorders.all(:disorder_type.not =>["DisEMBL Hotloops", "DisEMBL Coils"]).map{|k| k.disorder_id}
    dis_hl = self.disorders.first(:disorder_type => "DisEMBL Hotloops")
    dis_c = self.disorders.first(:disorder_type => "DisEMBL Coils")
    thread_array=[]
    thread_num.times do |i|
      thread_array[i] = Thread.new{
        while aa_array.length > 0 do
          aa = aa_array.pop
          # dis_sum = 0
          #           dvals = aa.disorder_values#DisorderValue.all(:aasequence_id => aa.id, :disorder_id=>dis_ids)
          #           dvs = dvals.map{|c| c.dvalue > 0.5 ? 1 : 0 }
          #           puts "DVals: " + dvs.join(',').to_s
          #           dis_sum = dvs.inject{|sum,x| sum + x }
          #           puts "Sum: " + dis_sum.to_s
          #           # if disorder_types.include?("DisEMBL Hotloops")
          #   dis_sum = dis_sum + 0.38 
          #   dis_num = dis_ids.length
          # end
          #end
          #if dvals.count > 0
          #  aa.disorder_consensus = dis_sum/dvals.count
          #else
          #  aa.disorder_consensus = 0
          #end
          aa.disorder_consensus = aa.disorder_values.all(:disorder_id => dis_ids, :dvalue.gte => 0.5 ).count#.avg(:dvalue)
          if !dis_hl.nil?
            if !aa.disorder_values.first(:disorder_id => dis_hl.id, :dvalue.gte => dis_hl.threshold).nil?
            aa.disorder_consensus  = aa.disorder_consensus + 1
            end
          end
          if !dis_c.nil?
            if !aa.disorder_values.first(:disorder_id => dis_c.id, :dvalue.gte => dis_c.threshold).nil?
            aa.disorder_consensus  = aa.disorder_consensus + 1
            end
          end
          aa.disorder_consensus = aa.disorder_consensus/self.disorders.count
          puts "Consensus: "+aa.disorder_consensus.to_s
          aa.save
        end
      }
    end
    thread_array.map{|t| t.join}
  end
  
  def self.calculate_disorder_consensus_for_types(ptype, disorder_types)
    Sequence.all(:seq_type=>ptype).each do |seq|
      seq.calculate_disorder_consensus(disorder_types)
    end
  end
  
  def generate_ronn
    unless self.disorders.first(:disorder_type=>"RONN")
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
    else
      puts "RONN Already Stored!***********************"
    end
   end

   def store_ronn(filepath)
     #create a new disorder object
     dis = Disorder.create(:seq_id => self.id, :disorder_type=>"RONN", :version=>1)
     self.disorders << dis
     self.save
     file = File.new(filepath, 'r')
     counter = 1
     aa_count = 0
     while (line = file.gets)
       #puts "#{counter}: #{line}"
       counter = counter + 1
       if counter > 2
         line_array = line.split
         if aa = AAsequence.first(:seq_id => self.id, :original_position=>aa_count)
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
           puts self.abrev_name + " Already Exists: #{seq.id}"
         else
           seq.generate_ronn
         end
       rescue
         puts seq.abrev_name + " Didn't Store: #{seq.id}"
       end
     end
   end

   def generate_disopred
     unless self.disorders.first(:disorder_type=>"Disopred")
       url ="http://bioinf.cs.ucl.ac.uk/psipred/submit"
       cur = Curl::Easy.new(url)
       post_params= "email=x58z545@montana.edu&program_disopred=1&sequence=#{self.sequence_one_line}&subject=#{self.abrev_name}_Disopred|#{self.id}&falseRate=2&output=opnone  "
       cur.http_post(post_params)
       puts post_params
       puts cur.body_str.to_s
       s =cur.body_str.to_s.split('<pre>')
       a = s[2].split("</pre>")
       filepath= "temp_data/#{self.abrev_name}_Disopred"
       f = File.new(filepath, "w+")
       f.write(a[0].to_s)
       f.close 
       puts filepath
       #self.store_ronn(filepath)
     else
       puts "Disopred Already Stored!***********************"
     end
    end

   def run_cornet
     unless !self.intra_residue_contacts.first.nil?
        puts "Submitting Cornet for: #{self.abrev_name}"
        url ="http://gpcr.biocomp.unibo.it/cgi/predictors/cornet/pred_cmapcgi.cgi"
        cur = Curl::Easy.new(url)
        post_params= "address=sean.b.cleveland@gmail.com&seqname=#{self.abrev_name}&db=SwissProt&text=#{self.sequence.gsub("\n","").gsub("\r","")}&filter=No"
        cur.http_post(post_params)
        puts post_params
         puts cur.body_str.to_s
     end
   end

   def run_svmseq
     unless !self.intra_residue_contacts.first.nil?
         puts "Submitting SVMSEQ for: #{self.abrev_name}"
         url ="http://zhanglab.ccmb.med.umich.edu/cgi-bin/SVMSEQ.pl"
         cur = Curl::Easy.new(url)
         #cur.multipart_form_post = true
         post_params= "REPLY-E-MAIL=sean.b.cleveland@gmail.com&TARGET-NAME=#{self.abrev_name}&SEQUENCE=#{self.sequence}"
         cur.http_post(post_params)
         puts post_params
         puts cur.body_str.to_s
     end
   end
   
   def run_dncon
     unless !self.intra_residue_contacts.first.nil?
          puts "Submitting DNCon for: #{self.abrev_name}"
          url ="http://iris.rnet.missouri.edu/cgi-bin/dncon/submit_job.cgi"
          cur = Curl::Easy.new(url)
          #cur.multipart_form_post = true
          post_params= "email_address=disiccapp@gmail.com&job_title=#{self.abrev_name}&protein_sequence=#{self.sequence}&to_take=l"
          cur.http_post(post_params)
          puts post_params
          puts cur.body_str.to_s
      end
   end


    

   def generate_pondr_fit
     unless self.disorders.first(:disorder_type=>"PONDR Fit")
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
     else
       puts "IUPred Already Stored!***********************"
    end
   end

   def store_pondr_fit(filepath)
     #create a new disorder object
     dis = Disorder.create(:seq_id => self.seq_id, :disorder_type=>"PONDR Fit", :version=>1)
     self.disorders << dis
     self.save
     file = File.new(filepath, 'r')
     counter = 1
     aa_count = 0
     while (line = file.gets)
       #puts "#{counter}: #{line}"
       counter = counter + 1
         line_array = line.split
         if aa = AAsequence.first(:seq_id => self.id, :original_position=>aa_count)
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
           puts seq.abrev_name + " Already Exists: #{seq.id}"
         else
           seq.generate_pondr_fit
         end
       rescue
         puts seq.abrev_name + " Didn't Store: #{seq.id}"
       end
     end
   end

   def generate_disembl
     if Disorder.first(:seq_id=>self.id, :disorder_type=>["DisEMBL Hotloops","DisEMBL Coils"])
        puts self.abrev_name + " DISEMBEL Already Exists: #{self.id}"
    else
       url = "http://dis.embl.de/cgiDict.py"
       cur = Curl::Easy.new(url)
       post_params= "key=process&SP_entry=&sequence_string=>#{self.sequence.gsub("\r",'').gsub("\n",'')}&smooth_frame=8&peak_frame=8&join_frame=4&fold_coils=1.20&fold_rem465=1.20&fold_hotloops=1.40&plot_title=&tango_PH=7.40&tango_T=298.15&tango_I=0.02&tango_TFE=0.00&fold_tango=1.00"
       puts post_params
       cur.http_post(post_params)
       s = cur.body_str.to_s.split('<th>Download predictions</th>')
       f = s[1].split(">smoothed scores</a>")
       a = f[0].to_s.gsub("\n<td><a href=",'').gsub('"','')
       file_url = "http://dis.embl.de/" + a
       threshold = cur.body_str.to_s.split("<th>Thresholds used</th>\n<td>")[1].split("</td>")[0];
       file_cur= Curl::Easy.new(file_url)
       file_cur.http_post()
       filepath= "temp_data/#{self.abrev_name}_Disembl"
       f = File.new(filepath, "w+")
       f.write(threshold + "\n")
       f.write(file_cur.body_str)
       f.close #this will be the smoothed scores file
       self.store_disembl(filepath)
     end
   end

   def store_disembl(filepath)
     #create a new disorder object
     dis_hl = Disorder.create(:seq_id => self.id, :disorder_type=>"DisEMBL Hotloops", :version=>1)
     self.disorders << dis_hl
     self.save
     dis_coil = Disorder.create(:seq_id => self.id, :disorder_type=>"DisEMBL Coils", :version=>1)
     self.disorders << dis_coil
     self.save
     file = File.new(filepath, 'r')
     counter = 1
     aa_count = 0
     while (line = file.gets)
       #puts "#{counter}: #{line}"
       counter = counter + 1
       if counter > 3
         line_array = line.split('      ')
         if aa = AAsequence.first(:seq_id => self.id, :original_position=>aa_count)
         #puts "Amino Acid -"+line_array[1]+ " : " + aa.amino_acid + " | " + aa_count.to_s
         #if aa.amino_acid == line_array[1]  
           DisorderValue.create(:disorder_id => dis_coil.disorder_id, :aasequence_id => aa.AAsequence_id, :dvalue=>line_array[0].to_f) 
           DisorderValue.create(:disorder_id => dis_hl.disorder_id, :aasequence_id => aa.AAsequence_id, :dvalue=>line_array[1].to_f) 
         end
         aa_count +=1
        elsif counter == 2
          coil = line.split("rem")[0].gsub("coils=",'').to_f
          loop = line.split("loops=")[1].to_f
          dis_hl.threshold = loop
          dis_coil.threshold = coil
          dis_hl.save
          dis_coil.save
        
       end
     
      
     end
   end

   def self.generate_disembl_for_all(seq_type)
     self.all(:seq_type => seq_type).each do |seq|
       begin
         if Disorder.first(:seq_id=>seq.seq_id, :disorder_type=>["DisEMBL Hotloops","DisEMBL Coils"])
           puts seq.abrev_name + " Already Exists: #{seq.id}"
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
  
  def import_conseq(alignment)
    seq = self      
    #open file that corresponds to this sequence
    puts filename = "temp_data/#{alignment.alignment_name}/conseq/#{seq.abrev_name}.conseq"
    if File.exists?(filename)
      seq.conseqs.destroy!
      puts "File exists"
      file = File.new(filename, "r")
      start_line = 13
      line_num = 1
      while (line = file.gets)
        if line_num > start_line
          results = line.split     
          puts results       
          cq = Conseq.new(:seq_id =>seq.seq_id,
                        :aasequence_id => seq.a_asequences.first(:original_position=>results[0].to_i-1).id,
                        :score => results[2].to_f,
                        :color => results[3].to_i,
                        :state => results[4],
                        :function => results[5],
                        :msa_data => results[6],
                        :residue_variety => results[7])
          cq.valid?
          puts cq.errors.inspect()
          cq.save
        end
        line_num +=1
      end #end while
    end #end if
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
  
  def generate_jalview_annotation_consensus()
    File.open("#{self.abrev_name}_#{self.seq_type}.gff", 'w') do |f1| 
      AAsequence.all(:seq_id => self.seq_id, :disorder_consensus.gte => 0.5).each do |aa|
        if aa.disorder_consensus > 0.5 && aa.disorder_consensus < 0.6
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
        f1.puts "None #{self.abrev_name} -1 #{aa.original_position} #{aa.original_position} #{feature_type}"
     end
   end
   #return jalview_string
  end
  
  def set_aasequence_defaults
      sql="UPDATE a_asequences SET disorder_consensus =0.0 WHERE seq_id=#{self.seq_id}"
      repository.adapter.execute(sql)
      sql="UPDATE a_asequences SET contact_consensus =0.0 WHERE seq_id=#{self.seq_id}"
      repository.adapter.execute(sql)
      sql="UPDATE a_asequences SET contact_positive_consensus =0.0 WHERE seq_id=#{self.seq_id}"
      repository.adapter.execute(sql)
  end
  
  def import_cornet(alignment)
    #open file that corresponds to this sequence
    seq = self
    puts filename = "temp_data/#{alignment.alignment_name}/cornet/#{seq.abrev_name}.cornet"
    if File.exists?(filename)
      seq.intra_residue_contacts(:type => "Cornet").destroy!
      puts "File exists"
      file = File.new(filename, "r")
      start_line = 17
      line_num = 1
      while (line = file.gets)
        if line_num > start_line
          results = line.split     
          puts results       
          intra = IntraResidueContact.new(:seq_id => seq.seq_id,
           :first_residue =>results[0].to_i,
           :second_residue =>results[1].to_i,
           :confidence => results[2].to_f,
           :type => "Cornet")
          intra.valid?
          puts intra.errors.inspect()
          intra.save
        end
        line_num +=1
      end #end while
    end #end if
  end
  
  def import_svmcon

      seq = self
      #open file that corresponds to this sequence
      puts filename = "#{self.abrev_name}_#{self.seq_type}.fasta.map"
      if File.exists?(filename)
        puts "File exists"
        file = File.new(filename, "r")
        start_line = 99999999999
        line_num = 1
        temp_lines = (self.a_asequences.count/50.0)
        lines = temp_lines.to_i + ( (temp_lines - temp_lines.to_i) > 0 ? 1 : 0)
        while (line = file.gets)
          if line.include?('Model')
            start_line = line_num + lines + 1
          end
          if line_num > start_line
            break if line == "\n"
            results = line.split
            puts "Result0:"+results[0]
            puts "Result1:"+results[1]
            puts results[0]
            position_one= results[0].to_i - 1
            puts results[1]
            position_two = results[1].to_i - 1
            d1 = results[2].to_i
            d2 = results[3].to_i
            confidence = results[4].to_f

            
            IntraResidueContact.create(:seq_id => self.id,
                        :first_residue => position_one,
                        :second_residue => position_two,
                        :d1 => mean_one,
                        :d2 => mean_two,
                        :confidence => correlation,
                        :seq_id => seq.seq_id,
                        :type => "SVMCon")
          end
          line_num +=1
        end #end while
      end #end if
  end
  
  #
  def evaluate_new_caps_distances(thread_num=200)
    caps_array = NewCap.all(:seq_id => self.seq_id).to_a
    thread_array=[]
    thread_num.times do |i|
      thread_array[i] = Thread.new{
        while caps_array.length > 0 do
          cap = caps_array.pop
          cap.eval_res_distances
        end
      }
    end
    thread_array.map{|t| t.join}
  end

  # new_caps_count
  # returns the number of caps records
  #
  def new_caps_count
    NewCap.all(:seq_id => self.seq_id).count
  end
  
  # new_caps_count_gt_twenty
  # returns the number of caps records that have pairs further than twenty amino acids away
  #
  def new_caps_count_gt_twenty
    NewCap.all(:seq_id => self.seq_id, :greater_than_twenty_away => true).count
  end


  # new_caps_count_gt_50
  # returns the number of caps records that have pairs further than 50 amino acids away
  #
  def new_caps_count_gt_50
    NewCap.all(:seq_id => self.seq_id, :greater_than_50_away => true).count
  end
  
  # new_caps_count_gt_100
  # returns the number of caps records that have pairs further than 1000 amino acids away
  #
  def new_caps_count_gt_100
    NewCap.all(:seq_id => self.seq_id, :greater_than_100_away => true).count
  end
  
  # new_caps_count_gt_100
  # returns the number of caps records that have pairs further than 1000 amino acids away
  #
  def new_caps_count_gt_200
    NewCap.all(:seq_id => self.seq_id, :greater_than_200_away => true).count
  end
  
  # new_caps_count_gt_100
  # returns the number of caps records that have pairs further than 1000 amino acids away
  #
  def new_caps_count_gt_300
    NewCap.all(:seq_id => self.seq_id, :greater_than_300_away => true).count
  end
  
  # new_caps_count_gt_100
  # returns the number of caps records that have pairs further than 1000 amino acids away
  #
  def new_caps_count_gt_400
    NewCap.all(:seq_id => self.seq_id, :greater_than_400_away => true).count
  end
  
  # new_caps_count_gt_100
  # returns the number of caps records that have pairs further than 1000 amino acids away
  #
  def new_caps_count_gt_500
    NewCap.all(:seq_id => self.seq_id, :greater_than_500_away => true).count
  end
  
  # new_caps_count_gt_100
  # returns the number of caps records that have pairs further than 1000 amino acids away
  #
  def new_caps_count_gt_600
    NewCap.all(:seq_id => self.seq_id, :greater_than_600_away => true).count
  end
  
  # new_caps_count_gt_100
  # returns the number of caps records that have pairs further than 1000 amino acids away
  #
  def new_caps_count_gt_700
    NewCap.all(:seq_id => self.seq_id, :greater_than_700_away => true).count
  end
  
  # new_caps_count_gt_100
  # returns the number of caps records that have pairs further than 1000 amino acids away
  #
  def new_caps_count_gt_800
    NewCap.all(:seq_id => self.seq_id, :greater_than_800_away => true).count
  end
  
  # new_caps_count_gt_100
  # returns the number of caps records that have pairs further than 1000 amino acids away
  #
  def new_caps_count_gt_900
    NewCap.all(:seq_id => self.seq_id, :greater_than_900_away => true).count
  end
  
  # new_caps_count_gt_100
  # returns the number of caps records that have pairs further than 1000 amino acids away
  #
  def new_caps_count_gt_1000
    NewCap.all(:seq_id => self.seq_id, :greater_than_1000_away => true).count
  end
  
  # new_caps_count_gt_100
  # returns the number of caps records that have pairs further than 1000 amino acids away
  #
  def new_caps_count_gt_1100
    NewCap.all(:seq_id => self.seq_id, :greater_than_1100_away => true).count
  end
  
  # new_caps_count_gt_100
  # returns the number of caps records that have pairs further than 1000 amino acids away
  #
  def new_caps_count_gt_1200
    NewCap.all(:seq_id => self.seq_id, :greater_than_1200_away => true).count
  end
  
  # new_caps_count_gt_100
  # returns the number of caps records that have pairs further than 1000 amino acids away
  #
  def new_caps_count_gt_1300
    NewCap.all(:seq_id => self.seq_id, :greater_than_1300_away => true).count
  end
  
  # new_caps_count_gt_100
  # returns the number of caps records that have pairs further than 1000 amino acids away
  #
  def new_caps_count_gt_1400
    NewCap.all(:seq_id => self.seq_id, :greater_than_1400_away => true).count
  end
  
  # new_caps_count_gt_100
  # returns the number of caps records that have pairs further than 1000 amino acids away
  #
  def new_caps_count_gt_1500
    NewCap.all(:seq_id => self.seq_id, :greater_than_1500_away => true).count
  end
  
  # new_caps_count_gt_100
  # returns the number of caps records that have pairs further than 1000 amino acids away
  #
  def new_caps_count_gt_1600
    NewCap.all(:seq_id => self.seq_id, :greater_than_1600_away => true).count
  end
  
  # new_caps_count_gt_100
  # returns the number of caps records that have pairs further than 1000 amino acids away
  #
  def new_caps_count_gt_1700
    NewCap.all(:seq_id => self.seq_id, :greater_than_1700_away => true).count
  end
  
  # new_caps_count_gt_100
  # returns the number of caps records that have pairs further than 1000 amino acids away
  #
  def new_caps_count_gt_1800
    NewCap.all(:seq_id => self.seq_id, :greater_than_1800_away => true).count
  end
  
  # new_caps_count_gt_100
  # returns the number of caps records that have pairs further than 1000 amino acids away
  #
  def new_caps_count_gt_1900
    NewCap.all(:seq_id => self.seq_id, :greater_than_1900_away => true).count
  end
  
  # new_caps_count_gt_100
  # returns the number of caps records that have pairs further than 1000 amino acids away
  #
  def new_caps_count_gt_2000
    NewCap.all(:seq_id => self.seq_id, :greater_than_2000_away => true).count
  end
  
  # new_caps_count_gt_100
  # returns the number of caps records that have pairs further than 1000 amino acids away
  #
  def new_caps_count_gt_2100
    NewCap.all(:seq_id => self.seq_id, :greater_than_2100_away => true).count
  end
  
  # new_caps_count_gt_100
  # returns the number of caps records that have pairs further than 1000 amino acids away
  #
  def new_caps_count_gt_2200
    NewCap.all(:seq_id => self.seq_id, :greater_than_2200_away => true).count
  end
  
  # new_caps_count_gt_100
  # returns the number of caps records that have pairs further than 1000 amino acids away
  #
  def new_caps_count_gt_2300
    NewCap.all(:seq_id => self.seq_id, :greater_than_2300_away => true).count
  end
  
  # new_caps_count_gt_100
  # returns the number of caps records that have pairs further than 1000 amino acids away
  #
  def new_caps_count_gt_2400
    NewCap.all(:seq_id => self.seq_id, :greater_than_2400_away => true).count
  end
  
  # new_caps_count_gt_100
  # returns the number of caps records that have pairs further than 1000 amino acids away
  #
  def new_caps_count_gt_2500
    NewCap.all(:seq_id => self.seq_id, :greater_than_2500_away => true).count
  end
end

