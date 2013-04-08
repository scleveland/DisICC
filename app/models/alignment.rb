class Alignment
  include DataMapper::Resource
  
  property :align_id, Serial, :key=>true
  property :seq_id, Integer, :required => true
  property :alignment_name, String, :required => true
  property :align_order, Integer, :required => true
  property :alignment_sequence, Text, :required => true, :default => ""
  property :fasta_title, Text, :required => false, :default => ""
  property :deleted_at, ParanoidDateTime
  
  # has n, :disorder_values
  has 1, :sequence, 'Sequence', :child_key => [:seq_id], :parent_key => [:seq_id]
  
  #validates_uniqueness_of :alignment_name
  
  #runs the score app to get percent identities for each sequence
  def run_score
    filename = self.generate_fasta_alignment_file
    string = "./lib/score_mac #{filename} temp_data/#{self.alignment_name}_res.txt temp_data/#{self.alignment_name}_dif.txt temp_data/#{self.alignment_name}_alignments.txt"
    puts string
    if system(string)
      
    end
  end
  
  def run_inter_caps_mac(dir_name, alignment1, alignment2)
    self.run_align_assess
    Dir.mkdir("temp_data/#{dir_name}") unless File.directory?("temp_data/#{dir_name}")
    alignments = Alignment.all(:alignment_name => self.alignment_name)
    alignment1.generate_pid_fasta_file_for_inter("temp_data/#{dir_name}")
    alignment2.generate_pid_fasta_file_for_inter("temp_data/#{dir_name}")
    system "./lib/comp_apps/caps/caps_mac -F temp_data/#{dir_name} --inter"    
  end
  
  def run_inter_caps_mac_without_fasta(dir_name)
    system "./lib/comp_apps/caps/caps_mac -F temp_data/#{dir_name} --inter"    
  end
  
  def alignments_to_positions_threaded(thread_num=4)
    alignments = Alignment.all(:alignment_name => self.alignment_name, :order=>[:align_order])
    alignment_array = []
    alignments.each do |alignment|
      alignment_array << alignment
    end
    thread_array=[]
    thread_num.times do |i|
      thread_array[i] = Thread.new{
        while alignment_array.length > 0 do
          align = alignment_array.pop
          AlignmentPosition.all(:alignment_id => align.align_id).destroy!
          align.alignment_to_positions()
        end
     }
    end
    thread_array.map{|t| t.join}
  end
  
  def alignment_to_positions
    puts self.sequence.abrev_name
    counter = 0
    aa_counter = 0
    sequence= self.sequence
    AlignmentPosition.all(:alignment_id => self.align_id).destroy!
    self.alignment_sequence.each_char do |aa|
      if aa != "-"
        #puts "counter: "+ counter.to_s + "  |aa :" + aa
        aaseq = sequence.a_asequences.first(:original_position=> aa_counter)
        if !aaseq.nil?
          AlignmentPosition.create(:alignment_id => self.align_id,
                            :position => counter,
                            :aasequence_id => aaseq.AAsequence_id)
        end
        aa_counter +=1
      end
      counter +=1
    end
  end
  
  def length_report
    alignments = Alignment.all(:alignment_name => self.alignment_name, :order=>[:align_order])
    alignments.each do |a|
      puts a.sequence.abrev_name + ": " + a.alignment_sequence.length.to_s + "  |seq_length :" + a.sequence.a_asequences.count.to_s
    end
  end
  
  def sequence
    Sequence.get(self.seq_id)
  end
  
  def sequences
    alignments = Alignment.all(:alignment_name => self.alignment_name, :order=>[:align_order])
    #seq_ids = alignments.map{|a| a.seq_id}
    #Sequence.all(:seq_id=>seq_ids)
    seq_array = []
    alignments.each do |a|
      seq_array << a.sequence
    end
    seq_array
  end
  
  
  #runs the disorder apps for every sequence in the alignment
  #then calculates the disorder consequence
  # @params id
  def run_consensus_disorder
    self.sequences.each do |sequence|
        sequence.run_and_store_disorder()
        sequence.calculate_disorder_consensus_threaded(300)
    end
  end
  
  #runs the AlignAssess app to get percent identities for each sequence in the alignment
  def run_align_assess
    filename = self.generate_fasta_alignment_file_for_all
    string = "./lib/AlignAssess_wShorterID #{filename} P"
    self.sequences.each do |s|
      PercentIdentity.all(:seq1_id=> s.seq_id, :alignment_name => self.alignment_name).destroy!
    end
    seq_array = Array.new
    if system(string)
      seq_id_array = self.sequences.map{|s| s.seq_id}
      new_filename = filename + "_assess"
      f = File.new(new_filename, "r")
      flag = false
      read_row= 999999999
      cur_row = 0
      while (line = f.gets)
        if cur_row > read_row  && flag
          if line == "\n"
            flag =false
          else
            puts line
            seq_array << line.split("\t")
          end
        elsif line == "Pair-wise %ID over shorter sequence:\n"
          flag=true
          read_row = cur_row + 2
        end
        cur_row +=1
      end
      range = seq_array.length - 1
      #seq_array.each do |row|
      for row_num in 0..range
        for i in 1..range+1#(row_num)      
          if row_num != i-1    
              PercentIdentity.first_or_create(:seq1_id=>seq_id_array[row_num],
                                                    :seq2_id=>seq_id_array[i-1],
                                                    :alignment_name => self.alignment_name,
                                                    :percent_id=>seq_array[row_num][i])
              PercentIdentity.first_or_create(:seq1_id=>seq_id_array[i-1],
                                                    :seq2_id=>seq_id_array[row_num],
                                                    :alignment_name => self.alignment_name,
                                                    :percent_id=>seq_array[row_num][i])
          end
          # print "[#{row_num}:#{i-1}=>#{row[i]}],"
        end
        #print "\n"
      end
    end
  end
  def run_caps_mac
    #self.run_align_assess
    Dir.mkdir("temp_data/#{self.alignment_name}") unless File.directory?("temp_data/#{self.alignment_name}")
    Dir.mkdir("temp_data/#{self.alignment_name}/caps") unless File.directory?("temp_data/#{self.alignment_name}/caps")
    alignments = Alignment.all(:alignment_name => self.alignment_name)
    alignments.each do |alignment|
      alignment.generate_pid_fasta_file_for_inter("temp_data/#{self.alignment_name}/caps")
    end
    system "./lib/comp_apps/caps/caps_mac -F temp_data/#{self.alignment_name}/caps --intra"
  end
  
  def run_caps_mac_with_old_fasta
    #self.run_align_assess
    # Dir.mkdir("temp_data/#{self.alignment_name}") unless File.directory?("temp_data/#{self.alignment_name}")
    #     Dir.mkdir("temp_data/#{self.alignment_name}/caps") unless File.directory?("temp_data/#{self.alignment_name}/caps")
    #     alignments = Alignment.all(:alignment_name => self.alignment_name)
    #     alignments.each do |alignment|
    #       alignment.generate_pid_fasta_file("temp_data/#{self.alignment_name}/caps")
    #     end
    system "./lib/comp_apps/caps/caps_mac -F temp_data/#{self.alignment_name}/caps --intra"
  end
  
  def run_caps
    self.run_align_assess
    Dir.mkdir("temp_data/#{self.alignment_name}") unless File.directory?("temp_data/#{self.alignment_name}")
    alignments = Alignment.all(:alignment_name => self.alignment_name)
    alignments.each do |alignment|
      alignment.generate_pid_fasta_file("temp_data/#{self.alignment_name}")
    end
    system "./lib/comp_apps/caps/caps -F temp_data/#{self.alignment_name} --intra"
  end
  
  def run_conseq(run_align_assess=true)
    if run_align_assess
      self.run_align_assess
    end
     Dir.mkdir("temp_data/#{self.alignment_name}") unless File.directory?("temp_data/#{self.alignment_name}")
     alignments = Alignment.all(:alignment_name => self.alignment_name)
     longest_alignment_length=0
     alignments.each do |alignment|
       if alignment.alignment_sequence.length > longest_alignment_length
          longest_alignment_length = alignment.alignment_sequence.length
       end
     end
     
     alignments.each do |alignment|
       if PercentIdentity.all(:seq1_id => alignment.sequence.id, :percent_id.gte => 19,:percent_id.lt => 90, :alignment_name => alignment.alignment_name).count > 9
         Dir.mkdir("temp_data/#{self.alignment_name}") unless File.directory?("temp_data/#{self.alignment_name}")
         alignment.generate_pid_fasta_file_for_inter("temp_data/#{self.alignment_name}", longest_alignment_length)
         #unless !self.conseq.first.nil?
           puts "Submitting Conseq for: #{alignment.sequence.abrev_name}"
           url ="http://conseq.tau.ac.il/cgi-bin/conseq_qsub.cgi"
           cur = Curl::Easy.new(url)
           cur.multipart_form_post = true
           post_params= "ALGORITHM=Bayes&email_address=sean.b.cleveland@gmail.com&msa_SEQNAME=#{alignment.sequence.abrev_name}&send_user_mail=yes"
           #cur.http_post(Curl::PostField.file('msa_FILE','temp_data/#{self.alignment_name}/#{self.alignment_name}_#{alignment.sequence.abrev_name}_pid.fasta'),post_params)
           myfile = File.open("temp_data/#{self.alignment_name}/#{self.alignment_name}_#{alignment.sequence.abrev_name}_pid.fasta", 'rb')
           # RestClient.post(url, :msa_FILE=>myfile, :ALGORITHM=>"Bayes",:email_address=>"sean.b.cleveland@gmail.com", :msa_SEQNAME=>alignment.sequence.abrev_name, :send_user_mail=>"yes"){ |response, request, result, &block|
           #   if [301, 302, 307].include? response.code
           #     myfile.close
           #     response.follow_redirection(request, result, &block)
           #   else
           #     myfile.close
           #     response.return!(request, result, &block)
           #   end
           # }
           system("curl -i -F ALGORITHM=Bayes -F email_address=sean.b.cleveland@gmail.com -F msa_SEQNAME=#{alignment.sequence.abrev_name} -F send_user_mail=yes -F msa_FILE=@temp_data/#{self.alignment_name}/#{self.alignment_name}_#{alignment.sequence.abrev_name}_pid.fasta http://conseq.tau.ac.il/cgi-bin/conseq_qsub.cgi")
           #cur.http_post(post_params)
           #puts post_params
           #puts cur.body_str.to_s
           #puts response.to_str
           break
         #end
       end
     end
  end
  
  def run_svmcon_threaded(thread_num = 4, create_fasta=true)
    alignments = Alignment.all(:alignment_name => self.alignment_name)
     alignment_array = []
     alignments.each do |alignment|
       alignment_array << alignment
     end
     thread_array=[]
     thread_num.times do |i|
       thread_array[i] = Thread.new{
         while alignment_array.length > 0 do
           alignment= alignment_array.pop
           if PercentIdentity.all(:seq1_id => alignment.sequence.seq_id, :percent_id.gte => 19,:percent_id.lt => 90, :alignment_name => alignment.alignment_name).count >= 10
             alignment.sequence.run_svmcon
           end
         end
      }
    end
    thread_array.map{|t| t.join}
   end
  
  
   def run_dncon_threaded(thread_num = 1)
     alignments = Alignment.all(:alignment_name => self.alignment_name)
      alignment_array = []
      alignments.each do |alignment|
        if PercentIdentity.all(:seq1_id => alignment.sequence.seq_id, :percent_id.gte => 19,:percent_id.lt => 90, :alignment_name => alignment.alignment_name).count >= 10
          alignment_array << alignment
        end
      end
      thread_array=[]
      thread_num.times do |i|
        thread_array[i] = Thread.new{
          while alignment_array.length > 0 do
            alignment= alignment_array.pop
            puts alignment.sequence.abrev_name + ": DNCON"
            alignment.sequence.run_dncon
          end
       }
     end
     thread_array.map{|t| t.join}
    end
  
   def run_nncon_threaded(thread_num = 4)
     alignments = Alignment.all(:alignment_name => self.alignment_name)
      alignment_array = []
      alignments.each do |alignment|
        alignment_array << alignment
      end
      Dir.mkdir("temp_data/#{self.alignment_name}/nncon") unless File.directory?("temp_data/#{self.alignment_name}/nncon")
      thread_array=[]
      thread_num.times do |i|
        thread_array[i] = Thread.new{
          while alignment_array.length > 0 do
            alignment= alignment_array.pop
            alignment.sequence.run_nncon("temp_data/#{self.alignment_name}/nncon")
          end
       }
     end
     thread_array.map{|t| t.join}
    end
    
    def run_cornet
     alignments = Alignment.all(:alignment_name => self.alignment_name)
     alignments.each do |alignment|
       alignment.sequence.run_cornet
     end
    end
    
  def run_caps_threaded(thread_num = 6)
     self.run_align_assess
     Dir.mkdir("temp_data/#{self.alignment_name}") unless File.directory?("temp_data/#{self.alignment_name}")
     alignments = Alignment.all(:alignment_name => self.alignment_name)
     longest_alignment_length=0
     alignments.each do |alignment|
       if alignment.alignment_sequence.length > longest_alignment_length
          longest_alignment_length = alignment.alignment_sequence.length
       end
     end
     alignment_array = []
     alignments.each do |alignment|
       if PercentIdentity.all(:seq1_id => alignment.sequence.id, :percent_id.gt => 25,:percent_id.lt => 90, :alignment_name => alignment.alignment_name).count > 9
         alignment_array << alignment.seq_id
         Dir.mkdir("temp_data/#{self.alignment_name}_#{alignment.seq_id}") unless File.directory?("temp_data/#{self.alignment_name}_#{alignment.seq_id}")
         alignment.generate_pid_fasta_file("temp_data/#{self.alignment_name}_#{alignment.seq_id}", longest_alignment_length)
       end
     end
     thread_array=[]
     thread_num.times do |i|
       thread_array[i] = Thread.new{
         while alignment_array.length > 0 do
           this_seq_id = alignment_array.pop
           puts "Starting Caps2 #{this_seq_id}"
           system "./lib/comp_apps/caps/caps -F temp_data/#{self.alignment_name}_#{this_seq_id} --intra"
           puts "Finsihed Caps2 for #{this_seq_id}"
         end
      }
    end
    thread_array.map{|t| t.join}
   end
  
   def run_caps_without_fasta_threaded(thread_num = 6)
      self.run_align_assess
      Dir.mkdir("temp_data/#{self.alignment_name}") unless File.directory?("temp_data/#{self.alignment_name}")
      alignments = Alignment.all(:alignment_name => self.alignment_name)
      alignment_array = []
      alignments.each do |alignment|
        if PercentIdentity.all(:seq1_id => alignment.sequence.id, :percent_id.gt => 25,:percent_id.lt => 90, :alignment_name => alignment.alignment_name).count > 9
          alignment_array << alignment.seq_id
        end
      end
      thread_array=[]
      thread_num.times do |i|
        thread_array[i] = Thread.new{
          while alignment_array.length > 0 do
            this_seq_id = alignment_array.pop
            puts "Starting Caps2 #{this_seq_id}"
            system "./lib/comp_apps/caps/caps -F temp_data/#{self.alignment_name}_#{this_seq_id} --intra"
            puts "Finsihed Caps2 for #{this_seq_id}"
          end
       }
     end
     thread_array.map{|t| t.join}
    end
  
    def run_caps_without_fasta
       #self.run_align_assess
       #Dir.mkdir("temp_data/#{self.alignment_name}") unless File.directory?("temp_data/#{self.alignment_name}")
       alignments = Alignment.all(:alignment_name => self.alignment_name)
       alignment_array = []
       alignments.each do |alignment|
         if PercentIdentity.all(:seq1_id => alignment.sequence.id, :percent_id.gt => 25,:percent_id.lt => 90, :alignment_name => alignment.alignment_name).count > 9
           alignment_array << alignment.seq_id
         end
       end
       #thread_array=[]
       #thread_num.times do |i|
       #   thread_array[i] = Thread.new{
           while alignment_array.length > 0 do
             this_seq_id = alignment_array.pop
             puts "Starting Caps2 #{this_seq_id}"
             system "./lib/comp_apps/caps2/caps2 -F temp_data/#{self.alignment_name}_#{this_seq_id} --intra"
             puts "Finsihed Caps2 for #{this_seq_id}"
           end
      #  }
      #end
      #thread_array.map{|t| t.join}
     end


  def self.run_inter_perl_caps(a1,a2,dir)
    Dir.mkdir("temp_data/#{dir}") unless File.directory?("temp_data/#{dir}")
    if PercentIdentity.all(:seq1_id => a1.sequence.id, :percent_id.gt => 19,:percent_id.lt => 90, :alignment_name => a1.alignment_name).count > 9
      a1.generate_pid_fasta_file_for_inter("temp_data/#{dir}")
    else
      puts "Not enough PIDs for #{a1.sequence.abrev_name}: #{a1.alignment_name}"
      return
    end
    if PercentIdentity.all(:seq1_id => a2.sequence.id, :percent_id.gt => 19,:percent_id.lt => 90, :alignment_name => a2.alignment_name).count > 9
      a2.generate_pid_fasta_file_for_inter("temp_data/#{dir}")
    else
      puts "Not enough PIDs for #{a2.sequence.abrev_name}: #{a2.alignment_name}"
      return
    end
    filepath = "temp_data/#{dir}/#{a1.alignment_name}_#{a1.alignment_name}_pid.ctl"
    f = File.new(filepath, "w+")
    ctl_string = "Input file1: temp_data/#{dir}/#{a1.alignment_name}_#{a1.sequence.abrev_name}_pid.fasta	* File containing sequence alignment for the first protein
    Input file2: temp_data/#{dir}/#{a2.alignment_name}_#{a2.sequence.abrev_name}_pid.fasta	* File containing sequence alignment for the second protein
    Out file1: temp_data/#{dir}/#{a1.alignment_name}_#{a2.alignment_name}_#{a2.sequence.abrev_name}_pid.out	* File where the output information should be stored
    Co-evolution analysis: 1			* (0) Intra-molecular; (1) Inter-protein
    Type of data 1: 0				* (0) amino acid alignment; (1) codon-based alignment
    Type of data 2: 0				* (0) amino acid alignment; (1) codon-based alignment
    3D test: 1						* Only applicable for intra-protein analysis: (0) perform test; (1) Test is not applicable
    Reference sequence 1: 1	* the order in the alignment of the sequences giving the real positions in the protein for 3D analyses
    Reference sequence 2: 1	* the order in the alignment of the sequences giving the real positions in the protein for 3D analyses
    3D file: 				* name of the file containing 3D coordinates
    Atom interval:1-1 		* Amino acid atoms for which 3D coordinates are available
    Significance test: 1				* (0) use threshold correlation; (1) random sampling
    Threshold R-value: 0.1			* Threshold value for the correlation coeficient
    Time correction: 1				* (0) no time correction; (1) weight correlations by the divergence time between sequences
    Time estimation: 0				* (0) use synonymous distances by Li 1993; (1) use Poisson-corrected amino acid distances
    Threshold alpha-value: 0.05		* This option valid only in case of random sampling
    Random sampling: 10000			* Use in case significance test option is 1
    Gaps: 2						* Remove all columns with gaps (0); Remove columns with a specified number of gaps (1); Do not remove columns with gaps(2)
    Minimum R: 0.1					* Minimum value of correlation coeficient to be considered for filtering
    GrSize: 3						* Maximum number of sites in the group permitted (given in percentage of protein length)"
    f.write(ctl_string)
    f.close
    puts "Starting Caps 1#{a1.sequence}"
    system "./lib/comp_apps/caps-perl/caps.pl #{filepath}"
    puts "Finsihed Caps 1 for #{a1.sequence}"
  end
  
  def run_aacon
    filename = self.generate_fasta_alignment_file_for_all
    Dir.mkdir("temp_data/#{self.alignment_name}") unless File.directory?("temp_data/#{self.alignment_name}")
    output_file = "temp_data/#{self.alignment_name}/conservation.txt"
    system "java -jar lib/conservation/compbio-conservation-1.0.jar  -i=#{filename} -m=KABAT  -o=#{output_file}"
    self.store_aacon
  end
  
  def store_aacon
    output_file = "temp_data/#{self.alignment_name}/conservation.txt"
    if File.exists?(output_file)
      puts "File exists**************************************************************************************************************"
      file = File.new(output_file, "r")
      line = file.gets
      Conservation.first_or_create(:alignment_id => self.align_id, :alignment_name => self.alignment_name, :conservation_results => line)
      file.close()
    else
     # self.run_aacon
    end
  end
  
  def self.run_inter_caps(a1,a2,dir)
    Dir.mkdir("temp_data/#{dir}") unless File.directory?("temp_data/#{dir}")
    if PercentIdentity.all(:seq1_id => a1.sequence.id, :percent_id.gt => 19,:percent_id.lt => 90, :alignment_name => a1.alignment_name).count > 9
      a1.generate_pid_fasta_file_for_inter("temp_data/#{dir}")
    else
      puts "Not enough PIDs for #{a1.sequence.abrev_name}: #{a1.alignment_name}"
      return
    end
    if PercentIdentity.all(:seq1_id => a2.sequence.id, :percent_id.gt => 19,:percent_id.lt => 90, :alignment_name => a2.alignment_name).count > 9
      a2.generate_pid_fasta_file_for_inter("temp_data/#{dir}")
    else
      puts "Not enough PIDs for #{a2.sequence.abrev_name}: #{a2.alignment_name}"
      return
    end
    puts "Starting Caps 1#{a1.sequence}"
    system "./lib/comp_apps/caps/caps_mac -F temp_data/#{dir} --inter"
    puts "Finsihed Caps 1 for #{a1.sequence}"
  end

     
  def run_perl_caps
    #self.run_align_assess
     #Dir.mkdir("temp_data/#{self.alignment_name}") unless File.directory?("temp_data/#{self.alignment_name}")
     thread_num = 2
     alignments = Alignment.all(:alignment_name => self.alignment_name)
     alignment_array = []
     alignments.each do |alignment|
       if PercentIdentity.all(:seq1_id => alignment.sequence.id, :percent_id.gt => 25,:percent_id.lt => 90, :alignment_name => alignment.alignment_name).count > 9
         alignment_array << alignment.seq_id
         Dir.mkdir("temp_data/#{self.alignment_name}") unless File.directory?("temp_data/#{self.alignment_name}")
         alignment.generate_pid_fasta_file("temp_data/#{self.alignment_name}")
       end
     end
    
     thread_array=[]
     thread_num.times do |i|
      thread_array[i] = Thread.new{
         while alignment_array.length > 0 do
           this_seq_id = alignment_array.pop
            filepath = "temp_data/#{self.alignment_name}/#{self.alignment_name}_#{Sequence.get(this_seq_id).abrev_name}_pid.ctl"
            f = File.new(filepath, "w+")
            ctl_string = "Input file1: temp_data/#{self.alignment_name}/#{self.alignment_name}_#{Sequence.get(this_seq_id).abrev_name}_pid.fasta	* File containing sequence alignment for the first protein
            Input file2: temp_data/#{self.alignment_name}/#{self.alignment_name}_#{Sequence.get(this_seq_id).abrev_name}_pid.fasta	* File containing sequence alignment for the second protein
            Out file1: temp_data/#{self.alignment_name}/#{self.alignment_name}_#{Sequence.get(this_seq_id).abrev_name}_pid.out	* File where the output information should be stored
            Co-evolution analysis: 0			* (0) Intra-molecular; (1) Inter-protein
            Type of data 1: 0				* (0) amino acid alignment; (1) codon-based alignment
            Type of data 2: 0				* (0) amino acid alignment; (1) codon-based alignment
            3D test: 1						* Only applicable for intra-protein analysis: (0) perform test; (1) Test is not applicable
            Reference sequence 1: 1	* the order in the alignment of the sequences giving the real positions in the protein for 3D analyses
            Reference sequence 2: 1	* the order in the alignment of the sequences giving the real positions in the protein for 3D analyses
            3D file: 				* name of the file containing 3D coordinates
            Atom interval:1-1 		* Amino acid atoms for which 3D coordinates are available
            Significance test: 1				* (0) use threshold correlation; (1) random sampling
            Threshold R-value: 0.1			* Threshold value for the correlation coeficient
            Time correction: 1				* (0) no time correction; (1) weight correlations by the divergence time between sequences
            Time estimation: 0				* (0) use synonymous distances by Li 1993; (1) use Poisson-corrected amino acid distances
            Threshold alpha-value: 0.05		* This option valid only in case of random sampling
            Random sampling: 10000			* Use in case significance test option is 1
            Gaps: 2						* Remove all columns with gaps (0); Remove columns with a specified number of gaps (1); Do not remove columns with gaps(2)
            Minimum R: 0.1					* Minimum value of correlation coeficient to be considered for filtering
            GrSize: 3						* Maximum number of sites in the group permitted (given in percentage of protein length)"
            f.write(ctl_string)
            f.close
           puts "Starting Caps 1#{this_seq_id}"
           system "./lib/comp_apps/caps-perl/caps.pl #{filepath}"
           puts "Finsihed Caps 1 for #{this_seq_id}"
         end
      }
    end
    thread_array.map{|t| t.join}
  end
  
  def run_perl_caps_without_fasta(thread_num=2)
    #self.run_align_assess
     #Dir.mkdir("temp_data/#{self.alignment_name}") unless File.directory?("temp_data/#{self.alignment_name}")
     alignments = Alignment.all(:alignment_name => self.alignment_name)
     alignment_array = []
     alignments.each do |alignment|
       if PercentIdentity.all(:seq1_id => alignment.sequence.id, :percent_id.gt => 25,:percent_id.lt => 90, :alignment_name => alignment.alignment_name).count > 9
         alignment_array << alignment.seq_id
       end
     end
    
     thread_array=[]
     thread_num.times do |i|
      thread_array[i] = Thread.new{
         while alignment_array.length > 0 do
           this_seq_id = alignment_array.pop
            filepath = "temp_data/#{self.alignment_name}/#{self.alignment_name}_#{Sequence.get(this_seq_id).abrev_name}_pid.ctl"
            f = File.new(filepath, "w+")
            ctl_string = "Input file1: temp_data/#{self.alignment_name}/#{self.alignment_name}_#{Sequence.get(this_seq_id).abrev_name}_pid.fasta	* File containing sequence alignment for the first protein
            Input file2: temp_data/#{self.alignment_name}/#{self.alignment_name}_#{Sequence.get(this_seq_id).abrev_name}_pid.fasta	* File containing sequence alignment for the second protein
            Out file1: temp_data/#{self.alignment_name}/#{self.alignment_name}_#{Sequence.get(this_seq_id).abrev_name}_pid.out	* File where the output information should be stored
            Co-evolution analysis: 0			* (0) Intra-molecular; (1) Inter-protein
            Type of data 1: 0				* (0) amino acid alignment; (1) codon-based alignment
            Type of data 2: 0				* (0) amino acid alignment; (1) codon-based alignment
            3D test: 1						* Only applicable for intra-protein analysis: (0) perform test; (1) Test is not applicable
            Reference sequence 1: 1	* the order in the alignment of the sequences giving the real positions in the protein for 3D analyses
            Reference sequence 2: 1	* the order in the alignment of the sequences giving the real positions in the protein for 3D analyses
            3D file: 				* name of the file containing 3D coordinates
            Atom interval:1-1 		* Amino acid atoms for which 3D coordinates are available
            Significance test: 1				* (0) use threshold correlation; (1) random sampling
            Threshold R-value: 0.1			* Threshold value for the correlation coeficient
            Time correction: 1				* (0) no time correction; (1) weight correlations by the divergence time between sequences
            Time estimation: 0				* (0) use synonymous distances by Li 1993; (1) use Poisson-corrected amino acid distances
            Threshold alpha-value: 0.05		* This option valid only in case of random sampling
            Random sampling: 10000			* Use in case significance test option is 1
            Gaps: 2						* Remove all columns with gaps (0); Remove columns with a specified number of gaps (1); Do not remove columns with gaps(2)
            Minimum R: 0.1					* Minimum value of correlation coeficient to be considered for filtering
            GrSize: 3						* Maximum number of sites in the group permitted (given in percentage of protein length)"
            f.write(ctl_string)
            f.close
           puts "Starting Caps 1#{this_seq_id}"
           system "./lib/comp_apps/caps-perl/caps.pl #{filepath}"
           puts "Finsihed Caps 1 for #{this_seq_id}"
         end
      }
    end
    thread_array.map{|t| t.join}
  end
  

  def import_rate4site(thread_num=4)
    seq_array = []
    self.sequences.each do |seq|
      seq_array << seq
    end
    thread_array=[]
    thread_num.times do |i|
    thread_array[i] = Thread.new{
       while seq_array.length > 0 do
          seq = seq_array.pop
          #puts filename = "temp_data/#{self.alignment_name}/#{self.alignment_name}_#{seq.abrev_name}_pid.fasta_xdet"#fasta.out"
          puts filename = "temp_data/#{self.alignment_name}/conseq/#{seq.abrev_name}.conseq"
          if File.exists?(filename)
            puts "File exists"
            file = File.new(filename, "r")
            while (line = file.gets)
                break if line == "\n"
                results = line.split
                puts "Postion:"+results[0]+ "| Conservation:" + results[4] + "| Correlation:" + results[8]
                xd = Conseq.new(
                  :aasequence_id =>AlignmentPosition.first(:position=> results[0].to_i-1,:alignment_id=>self.align_id).aasequence_id,
                  :seq_id => seq.seq_id,
                  :score => results[2],
                  :color => results[2],
                  :state => 0,
                  :function => 0,
                  :msa_data => results[5],
                  :residue_variety => 0,
                )  
                xd.valid?
                puts xd.errors.inspect()
                xd.save
            end #end while
          end #end if
        end
      }
    end
    thread_array.map{|t| t.join}
  end
  
  def import_rate4site_single()
   seq = self.sequence
   #puts filename = "temp_data/#{self.alignment_name}/#{self.alignment_name}_#{seq.abrev_name}_pid.fasta_xdet"#fasta.out"
   puts filename = "temp_data/#{self.alignment_name}/#{seq.abrev_name}.conseq"
   if File.exists?(filename)
     puts "File exists"
     file = File.new(filename, "r")
     while (line = file.gets)
         break if line == "\n"
         results = line.split
         puts "Postion:"+results[0]+ "| Conservation:" + results[4] + "| Correlation:" + results[8]
         xd = Conseq.new(
           :aasequence_id =>AlignmentPosition.first(:position=> results[0].to_i-1,:alignment_id=>self.align_id).aasequence_id,
           :seq_id => seq.seq_id,
           :score => results[2],
           :color => results[2],
           :state => 0,
           :function => 0,
           :msa_data => results[5],
           :residue_variety => 0,
         )  
         xd.valid?
         puts xd.errors.inspect()
         xd.save
     end #end while
   end #end if
  end
  
  def import_xdet(thread_num=4)
    seq_array = []
    self.sequences.each do |seq|
      seq_array << seq
    end
    thread_array=[]
    thread_num.times do |i|
    thread_array[i] = Thread.new{
       while seq_array.length > 0 do
          seq = seq_array.pop
          puts filename = "temp_data/#{self.alignment_name}/#{self.alignment_name}_#{seq.abrev_name}_pid.fasta_xdet"#fasta.out"
          if File.exists?(filename)
            puts "File exists"
            puts "File exists"
            puts "File exists**************************************************************************************************************"
            file = File.new(filename, "r")
            while (line = file.gets)
                break if line == "\n"
                results = line.split
                puts "Postion:"+results[0]+ "| Conservation:" + results[4] + "| Correlation:" + results[8]
                begin
                xd = Xdet.new(
                  :aasequence_id =>AlignmentPosition.first(:position=> results[0].to_i-1,:alignment_id=>self.align_id).aasequence_id,
                  :conservation => results[8].to_f, 
                  :correlation => results[4].to_f, 
                  :seq_id => seq.seq_id
                )  
                xd.valid?
                puts xd.errors.inspect()
                xd.save
                rescue Exception => e
                  puts "Something went not right. It went wrong in fact. It went uncorrect.!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                  puts e.message
                end
            end #end while
          end #end if
        end
      }
    end
    thread_array.map{|t| t.join}
  end
  
  def import_xdet_for_seq(thread_num=4)
    seq_array = []
    seq_array << self.sequence
    thread_array=[]
    thread_num.times do |i|
    thread_array[i] = Thread.new{
       while seq_array.length > 0 do
          seq = seq_array.pop
          puts filename = "temp_data/#{self.alignment_name}/#{self.alignment_name}_#{seq.abrev_name}_pid.fasta_xdet"#fasta.out"
          if File.exists?(filename)
            puts "File exists"
            puts "File exists"
            puts "File exists**************************************************************************************************************"
            file = File.new(filename, "r")
            while (line = file.gets)
                break if line == "\n"
                results = line.split
                puts "Postion:"+results[0]+ "| Conservation:" + results[4] + "| Correlation:" + results[8]
                begin
                xd = Xdet.new(
                  :aasequence_id =>AlignmentPosition.first(:position=> results[0].to_i-1,:alignment_id=>self.align_id).aasequence_id,
                  :conservation => results[8].to_f, 
                  :correlation => results[4].to_f, 
                  :seq_id => seq.seq_id
                )  
                xd.valid?
                puts xd.errors.inspect()
                xd.save
                rescue Exception => e
                  puts "Something went not right. It went wrong in fact. It went uncorrect.!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                  puts e.message
                end
            end #end while
          end #end if
        end
      }
    end
    thread_array.map{|t| t.join}
  end

  def import_xdet_new(thread_num=4)
    alignments = Alignment.all(:alignment_name => self.alignment_name)
    alignment_array = []
    alignments.each do |a|
      alignment_array << a
    end
    thread_array=[]
    thread_num.times do |i|
    thread_array[i] = Thread.new{
       while alignment_array.length > 0 do
          alignment = alignment_array.pop
           if PercentIdentity.all(:seq1_id => alignment.sequence.id, :percent_id.gte => 19,:percent_id.lt => 90, :alignment_name => alignment.alignment_name).count > 9
            puts filename = "temp_data/#{alignment.alignment_name}/#{alignment.alignment_name}_#{alignment.sequence.abrev_name}_pid.fasta_xdet"#fasta.out" 
            if File.exists?(filename)
              puts "File exists"
              puts "File exists"
              puts "File exists**************************************************************************************************************"
              file = File.new(filename, "r")
              while (line = file.gets)
                  break if line == "\n"
                  results = line.split
                  puts "Postion:"+results[0]+ "| Conservation:" + results[4] + "| Correlation:" + results[8]
                  begin
                  xd = Xdet.new(
                    :aasequence_id =>AlignmentPosition.first(:position=> results[0].to_i-1,:alignment_id=>alignment.align_id).aasequence_id,
                    #:position => results[0].to_i-1,
                    :conservation => results[8].to_f, 
                    :correlation => results[4].to_f, 
                    :seq_id => alignment.sequence.seq_id
                  )  
                  xd.valid?
                  puts xd.errors.inspect()
                  xd.save
                  rescue Exception => e
                    puts "Something went not right. It went wrong in fact. It went uncorrect.!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                    puts e.message
                  end
              end #end while
            end #end if
          end
        end
      }
    end
    thread_array.map{|t| t.join}
  end
  
  def import_xdet_missing(thread_num=4)
    alignments = Alignment.all(:alignment_name => self.alignment_name)
    alignment_array = []
    alignments.each do |a|
      if Xdet.all(:seq_id => a.seq_id).count == 0 && PercentIdentity.all(:seq1_id => a.seq_id, :percent_id.gte => 19, :order =>[:percent_id.desc],:alignment_name=>a.alignment_name).count > 9
        alignment_array << a
      end
    end
    #thread_array=[]
    #thread_num.times do |i|
    #thread_array[i] = Thread.new{
       while alignment_array.length > 0 do
          alignment = alignment_array.pop
          puts filename = "temp_data/#{alignment.alignment_name}/#{alignment.alignment_name}_#{alignment.sequence.abrev_name}_pid.fasta_xdet"#fasta.out"
          if File.exists?(filename)
            puts "File exists"
            puts "File exists"
            puts "File exists**************************************************************************************************************"
            file = File.new(filename, "r")
            while (line = file.gets)
                break if line == "\n"
                results = line.split
                puts "Postion:"+results[0]+ "| Conservation:" + results[4] + "| Correlation:" + results[8]
                begin
                xd = Xdet.new(
                  :aasequence_id =>AlignmentPosition.first(:position=> results[0].to_i-1,:alignment_id=>alignment.align_id).aasequence_id,
                  :conservation => results[8].to_f, 
                  :correlation => results[4].to_f, 
                  :seq_id => alignment.sequence.seq_id
                )  
                xd.valid?
                puts xd.errors.inspect()
                xd.save
                rescue Exception => e
                  puts "Something went not right. It went wrong in fact. It went uncorrect.!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                  puts e.message
                end
            end #end while
          end #end if
        end
    #  }
    #end
    #thread_array.map{|t| t.join}
  end

  def run_xdet
    #self.run_align_assess
    Dir.mkdir("#{dir}") unless File.directory?("#{dir}")
    Dir.mkdir("temp_data/#{self.alignment_name}") unless File.directory?("temp_data/#{self.alignment_name}")
    alignments = Alignment.all(:alignment_name => self.alignment_name)
    alignments.each do |alignment|
      filename= alignment.generate_pid_fasta_file_for_inter("temp_data/#{self.alignment_name}")
      system "./lib/comp_apps/XDet/xdet_linux32 #{filename} ~/Rails/DisICC/lib/comp_apps/XDet/Maxhom_McLachlan.metric >> #{filename}_xdet"
    end
  end
  
  def run_xdet_threaded(dir="temp_data", thread_num=4)
    #self.run_align_assess
    Dir.mkdir("#{dir}") unless File.directory?("#{dir}")
    Dir.mkdir("#{dir}/#{self.alignment_name}") unless File.directory?("#{dir}/#{self.alignment_name}")
    alignments = Alignment.all(:alignment_name => self.alignment_name)
    alignment_array = []
    alignments.each do |a|
      alignment_array << a
    end
    thread_array=[]
     thread_num.times do |i|
       thread_array[i] = Thread.new{
         while alignment_array.length > 0 do
          alignment = alignment_array.pop
          #filename= alignment.generate_pid_fasta_file("temp_data/#{self.alignment_name}")
          filename= alignment.generate_pid_fasta_file_for_inter("#{dir}/#{self.alignment_name}")
          system "./lib/comp_apps/XDet/xdet_linux32 #{filename} ~/Rails/DisICC/lib/comp_apps/XDet/Maxhom_McLachlan.metric >> #{filename}_xdet"
        end
      }
    end
    thread_array.map{|t| t.join}
  end
  
  def run_xdet_without_fasta_threaded(thread_num=4)
    #self.run_align_assess
    Dir.mkdir("temp_data/#{self.alignment_name}") unless File.directory?("temp_data/#{self.alignment_name}")
    alignments = Alignment.all(:alignment_name => self.alignment_name)
    alignment_array = []
    alignments.each do |a|
      alignment_array << a
    end
    thread_array=[]
     thread_num.times do |i|
       thread_array[i] = Thread.new{
         while alignment_array.length > 0 do
          alignment = alignment_array.pop
          #filename= alignment.generate_pid_fasta_file("temp_data/#{self.alignment_name}")
          #filename= alignment.generate_pid_fasta_file_for_inter("temp_data/#{self.alignment_name}")
          filename = "temp_data/#{self.alignment_name}/#{self.alignment_name}_#{alignment.sequence.abrev_name}_pid.fasta"
          system "./lib/comp_apps/XDet/xdet_linux32 #{filename} ~/Rails/DisICC/lib/comp_apps/XDet/Maxhom_McLachlan.metric >> #{filename}_xdet"
        end
      }
    end
    thread_array.map{|t| t.join}
  end
  

  def run_rate4site_single
    Dir.mkdir("temp_data/#{self.alignment_name}") unless File.directory?("temp_data/#{self.alignment_name}")
    filename = "temp_data/#{self.alignment_name}/#{self.alignment_name}_#{self.sequence.abrev_name}_pid.fasta"
    system "./lib/comp_apps/conseq/rate4site_fast -s #{filename} -o #{self.sequence.abrev_name}.conseq"
  end

  def run_rate4site
    #self.run_align_assess
    Dir.mkdir("temp_data/#{self.alignment_name}") unless File.directory?("temp_data/#{self.alignment_name}")
    alignments = Alignment.all(:alignment_name => self.alignment_name)
    alignments.each do |alignment|
      #filename= alignment.generate_pid_fasta_file_for_inter("temp_data/#{self.alignment_name}")
      filename = "temp_data/#{self.alignment_name}/#{self.alignment_name}_#{alignment.sequence.abrev_name}_pid.fasta"
      system "./lib/comp_apps/conseq/rate4site_fast -s #{filename} -o temp_data/#{self.alignment_name}/#{alignment.sequence.abrev_name}.conseq"
    end
  end
  
  def generate_pid_fasta_files
    self.run_align_assess
    self.sequences.each do |seq|
      fasta_string=""
      pids = PercentIdentity.all(:seq1_id => seq.seq_id, :percent_id.gte => 19, :order =>[:percent_id.desc])
      if pids.count > 10
        fasta_string= Alignment.first(:alignment_name => self.alignment_name, :seq_id=>seq.seq_id).fasta_alignment_string
        puts seq.abrev_name+":"+pids.count.to_s
        pids.each do |pid|
          if pid.seq2_id != seq.seq_id
            print Sequence.get(pid.seq2_id).abrev_name + ":" + pid.percent_id.to_s + ","
            fasta_string = fasta_string + Alignment.first(:alignment_name=>pid.alignment_name, :seq_id=>pid.seq2_id).fasta_alignment_string("pid:#{pid.percent_id}")
          end
        end
        filepath = "temp_data/"+self.alignment_name+"_"+seq.abrev_name+"_pid.fasta"
        f = File.new(filepath, "w+")
        f.write(fasta_string)
        f.close
        print "\n"
      end
    end
  end
 
  def generate_pid_inter_fasta_files(dir="temp_data")
    alignments = Alignment.all(:alignment_name => self.alignment_name)
    alignments.each do |a|
      if PercentIdentity.all(:seq1_id => a.sequence.seq_id, :percent_id.gte => 19, :order =>[:percent_id.desc], :alignment_name => self.alignment_name).count > 9
        a.generate_pid_fasta_file_for_inter(dir)
      end
    end
  end
  #assumes that run_align_asses has already occured
   def generate_pid_fasta_file_for_inter(dir="temp_data", longest_alignment_length=0)
     fasta_string=""
     Dir.mkdir("#{dir}") unless File.directory?("#{dir}")
     seq = Sequence.get(self.seq_id)
     pids = PercentIdentity.all(:seq1_id => self.seq_id, :percent_id.gte => 19, :order =>[:percent_id.desc],:unique=>true)
     fasta_string= Alignment.first(:alignment_name => self.alignment_name, :seq_id=>self.seq_id).fasta_alignment_string_inter()
     puts seq.abrev_name+":"+pids.count.to_s
     puts pids.map{|p| p.seq2_sequence.seq_name}.join(',')
     seq_array = pids.map{|p| p.seq2_id}.uniq  
     seq_array.each do |seq_id|
       if seq_id != seq.id
         if Alignment.first(:alignment_name=>self.alignment_name, :seq_id=>seq_id)
           print Sequence.get(seq_id).abrev_name + ":" + pids.first(:seq2_id =>seq_id).percent_id.to_s + ","
          fasta_string = fasta_string + Alignment.first(:alignment_name=>self.alignment_name, :seq_id=>seq_id).fasta_alignment_string_inter()
         end
       end
     end
     puts ""
     filepath = "#{dir}/"+self.alignment_name+"_"+seq.abrev_name+"_pid.fasta"
     f = File.new(filepath, "w+")
     f.write(fasta_string)
     f.close
     filepath
   end
 
  
  #assumes that run_align_asses has already occured
  def generate_pid_fasta_file(dir="temp_data", longest_alignment_length=0)
    fasta_string=""
    seq = Sequence.get(self.seq_id)
    pids = PercentIdentity.all(:seq1_id => self.seq_id, :percent_id.gte => 19, :order =>[:percent_id.desc],:unique=>true)
    fasta_string= Alignment.first(:alignment_name => self.alignment_name, :seq_id=>self.seq_id).fasta_alignment_string("",longest_alignment_length)
    puts seq.abrev_name+":"+pids.count.to_s
    puts pids.map{|p| p.seq2_sequence.seq_name}.join(',')
    pids.each do |pid|
      if pid.seq2_id != seq.id
        print Sequence.get(pid.seq2_id).abrev_name + ":" + pid.percent_id.to_s + ","
        fasta_string = fasta_string + Alignment.first(:alignment_name=>pid.alignment_name, :seq_id=>pid.seq2_id).fasta_alignment_string("pid:#{pid.percent_id}",longest_alignment_length)
      end
    end
    puts ""
    filepath = "#{dir}/"+self.alignment_name+"_"+seq.abrev_name+"_pid.fasta"
    f = File.new(filepath, "w+")
    f.write(fasta_string)
    f.close
    filepath
  end
  
  
  #assumes that run_align_asses has already occured
  def generate_pid_fasta_file_over(dir="temp_data", longest_alignment_length=0, over_num =9)
    fasta_string=""
    seq = Sequence.get(self.seq_id)
    pids = PercentIdentity.all(:seq1_id => self.seq_id, :percent_id.gte => 19, :order =>[:percent_id.desc],:unique=>true)
    if pids.count > over_num
      fasta_string= Alignment.first(:alignment_name => self.alignment_name, :seq_id=>self.seq_id).fasta_alignment_string("",longest_alignment_length)
      puts seq.abrev_name+":"+pids.count.to_s
      puts pids.map{|p| p.seq2_sequence.seq_name}.join(',')
      pids.each do |pid|
        if pid.seq2_id != seq.id
          print Sequence.get(pid.seq2_id).abrev_name + ":" + pid.percent_id.to_s + ","
          fasta_string = fasta_string + Alignment.first(:alignment_name=>pid.alignment_name, :seq_id=>pid.seq2_id).fasta_alignment_string()
        end
      end
      puts ""
      filepath = "#{dir}/"+self.alignment_name+"_"+seq.abrev_name+"_pid.fasta"
      f = File.new(filepath, "w+")
      f.write(fasta_string)
      f.close
      filepath
    end
  end
  
  def update_alignment_sequence
    cur_position = 0   
    fasta_string = ""
    AlignmentPosition.all(:alignment_id => self.align_id, :order => [:alignmnet_position_id.asc]).each do |position|
      if position.position == cur_position
         fasta_string = fasta_string + AAsequence.first(:AAsequence_id => position.aasequence_id).amino_acid
      else
         while position.position > cur_position
                      fasta_string = fasta_string +"-"
                      cur_position += 1
         end
        fasta_string = fasta_string + AAsequence.first(:AAsequence_id => position.aasequence_id).amino_acid
      end
      cur_position += 1         
    end
    self.alignment_seqence = fasta_string
    self.save
  end
  def fasta_alignment_string_inter(longest_alignment_length=0)
    if self.alignment_sequence.empty?
      self.update_alignment_sequence
    end
    extra_gaps = longest_alignment_length - self.alignment_sequence.length
    extra_str=""
    unless extra_gaps < 0
      extra_gaps.times do |i|
        extra_str = extra_str + "-"
      end
    end
    seq = Sequence.get(self.seq_id)
    fasta_string = ">"+seq.abrev_name + "\n" + self.alignment_sequence.gsub("\n","") + "\n"
  end
  
  def fasta_alignment_string(extra_string="",longest_alignment_length=0)
    if self.alignment_sequence.empty?
      self.update_alignment_sequence
    end
    extra_gaps = longest_alignment_length - self.alignment_sequence.length
    extra_str=""
    unless extra_gaps < 0
      extra_gaps.times do |i|
        extra_str = extra_str + "-"
      end
    end
    seq = Sequence.get(self.seq_id)
    fasta_string = ">"+seq.abrev_name+" | "+seq.seq_name+"|"+seq.seq_type+"|"+seq.seq_accession+"|"+extra_string+"\n" + self.alignment_sequence+ extra_str+"\n"
  end
  
  def generate_fasta_alignment_file_for_all(filename="", dir="temp_data")
    Dir.mkdir("#{dir}") unless File.directory?("#{dir}")
    alignments = Alignment.all(:alignment_name=> self.alignment_name, :order =>[:align_id])
    if filename.empty?
      filepath = "#{dir}/"+self.alignment_name+"_alignment"+Time.now.to_i.to_s+".fasta"
    else
      filepath = "#{dir}/#{filename}"
    end
    f = File.new(filepath, "w+")
    alignments.each do |alignment|
      f.write(alignment.fasta_alignment_string_inter)
    end
    f.close
    filepath
  end
  
  def import_svmcon

    self.sequences.each do |seq|
      seq.import_svmcon
    end #end sequences.each
  end

  def self.import_inter_caps(a1,a2,path)
  seq1 = a1.sequence
  seq2 = a2.sequence
  a1_first = Alignment.first(:alignment_name => a1.alignment_name)
  a2_first = Alignment.first(:alignment_name => a2.alignment_name)
    #open file that corresponds to this sequence
    puts filename = path#"temp_data/#{dir}/#{a1.alignment_name}_#{a2.alignment_name}_#{a2.sequence.abrev_name}_pid.out"
    if File.exists?(filename)
      puts "File exists"
      file = File.new(filename, "r")
      start_line = 99999999999
      line_num = 1
      while (line = file.gets)
        if line.include?('Posicion')
          start_line = line_num + 1
        end
        if line_num > start_line
          break if line == "\n"
          results = line.split
          puts "Result0:"+results[0]
          puts "Result1:"+results[1]
          puts results[0].split('(')[1]
          position_one= results[0].split('(')[1].gsub(')','').to_i - 1
          aasequence_one = AAsequence.first(:seq_id=> seq1.seq_id, :original_position=>position_one)
          puts results[1].split('(')[1]
          position_two = results[1].split('(')[1].gsub(')','').to_i - 1
          aasequence_two = AAsequence.first(:seq_id=> seq2.seq_id, :original_position=>position_two)
          mean_one = results[2].to_f
          mean_two = results[3].to_f
          correlation = results[4].to_f
          begin
          InterCap.create(:aasequence1_id => aasequence_one.AAsequence_id,
                      :aasequence2_id => aasequence_two.AAsequence_id,
                      :position_one => position_one,
                      :position_two => position_two,
                      :mean_one => mean_one,
                      :mean_two => mean_two,
                      :correlation => correlation,
                      :seq1_id => seq1.seq_id,
                      :seq2_id => seq2.seq_id,
                      :alignment1_id => a1_first.align_id,
                      :alignment2_id => a2_first.align_id )
          rescue Exception => e
           puts e
          end
        end
        line_num +=1
      end #end while
    end #end if
  end

  def import_caps
    #find correct directory
    #dir_name = "#{root_dir}/#{self.alignment_name}"
    #for each sequence in the alignment
    self.sequences.each do |seq|
      #open file that corresponds to this sequence
      puts filename = "temp_data/#{self.alignment_name}_#{seq.abrev_name}_pid.out"#fasta.out"
      if File.exists?(filename)
        puts "File exists"
        file = File.new(filename, "r")
        start_line = 99999999999
        line_num = 1
        while (line = file.gets)
          if line.include?('Posicion')
            start_line = line_num + 1
          end
          if line_num > start_line
            break if line == "\n"
            results = line.split
            puts "Result0:"+results[0]
            puts "Result1:"+results[1]
            puts results[0].split('(')[1]
            position_one= results[0].split('(')[1].gsub(')','').to_i - 1
            aasequence_one = AAsequence.first(:seq_id=> seq.seq_id, :original_position=>position_one)
            puts results[1].split('(')[1]
            position_two = results[1].split('(')[1].gsub(')','').to_i - 1
            aasequence_two = AAsequence.first(:seq_id=> seq.seq_id, :original_position=>position_two)
            mean_one = results[2].to_f
            mean_two = results[3].to_f
            correlation = results[4].to_f
            NewCap.create(:aasequence_id => aasequence_one.AAsequence_id,
                        :position_one => position_one,
                        :position_two => position_two,
                        :mean_one => mean_one,
                        :mean_two => mean_two,
                        :correlation => correlation,
                        :seq_id => seq.seq_id )
            # NewCap.create(:aasequence_id => aasequence_one.AAsequence_id,
            #             :position_one => position_two,
            #             :position_two => position_one,
            #             :mean_one => mean_one,
            #             :mean_two => mean_two,
            #             :correlation => correlation,
            #             :seq_id => seq.seq_id )
          end
          line_num +=1
        end #end while
      end #end if
    end #end sequences.each
  end

  def generate_aasequences(thread_num=65)
    seq_array=[]
    self.sequences.each do |seq|
      seq_array << seq
    end
    thread_array=[]
    thread_num.times do |i|
      thread_array[i] = Thread.new{
         while seq_array.length > 0 do
            seq = seq_array.pop
            puts seq.abrev_name
            seq.generate_aasequences
        end
      }
    end
    thread_array.map{|t| t.join}
  end
  
  def run_and_store_disorder(thread_num=65)
    seq_array=[]
    self.sequences.each do |seq|
      seq_array << seq
    end
    thread_array=[]
    thread_num.times do |i|
      thread_array[i] = Thread.new{
         while seq_array.length > 0 do
            seq = seq_array.pop
            puts seq.abrev_name
            seq.run_and_store_disorder
        end
      }
    end
    thread_array.map{|t| t.join}
  end
  
  def store_disorder(thread_num=65)
    seq_array=[]
    self.sequences.each do |seq|
      seq_array << seq
    end
    thread_array=[]
    thread_num.times do |i|
      thread_array[i] = Thread.new{
         while seq_array.length > 0 do
            seq = seq_array.pop
            puts seq.abrev_name
            seq.store_disorder
        end
      }
    end
    thread_array.map{|t| t.join}
  end
  
  def undelete_aasequences
    #AlignmentPosition.all(:alignment_id => self.align_id).each do |ap|
    #AlignmentPosition.all(:alignment_id => self.align_id).each do |apos|
      sql = "UPDATE a_asequences SET deleted_at = NULL WHERE a_asequence_id in( #{AlignmentPosition.all(:alignment_id => self.align_id).map{|ap| ap.aasequence_id}.join(',')})"
      repository.adapter.execute(sql)
    #end 
  end

  # def delete_disorders
  #   self.sequences.each do |seq|
  #     puts seq.abrev_name
  #     seq.disorders.destroy!
  #   end
  # end

  def delete_caps
    self.sequences.each do |seq|
      puts seq.abrev_name
      NewCap.all(:seq_id => seq.seq_id).destroy!
    end
  end
  
  def delete_conseq
    self.sequences.each do |seq|
      puts seq.abrev_name
     Conseq.all(:seq_id => seq.seq_id).destroy!
    end
  end

  def delete_xdet
    self.sequences.each do |seq|
      puts seq.abrev_name
      Xdet.all(:seq_id => seq.seq_id).destroy!
    end
  end


  def import_caps_threaded(dir = "temp_data/", thread_num=4)
    seq_array = []
    self.sequences.each do |seq|
      if NewCap.all(:seq_id => seq.seq_id).count <= 0 
         seq_array << seq
         puts seq.abrev_name
      else
        puts seq.seq_name + "| " + NewCap.all(:seq_id => seq.seq_id).count.to_s
      end
    end
    thread_array=[]
    thread_num.times do |i|
      thread_array[i] = Thread.new{
         while seq_array.length > 0 do
            seq = seq_array.pop
            if PercentIdentity.all(:seq1_id => seq.seq_id, :percent_id.gte => 19,:percent_id.lt => 90, :alignment_name => self.alignment_name).count > 9
              puts filename = "#{dir}#{self.alignment_name}_#{seq.abrev_name}_pid.fasta.out"
              #temp_data/#{self.alignment_name}/#{self.alignment_name}_#{Sequence.get(this_seq_id).abrev_name}_pid.fasta
              if File.exists?(filename)
                puts "File exists"
                file = File.new(filename, "r")
                start_line = 99999999999
                line_num = 1
                while (line = file.gets)
                  if line.include?('Posicion') ||  line.include?('Position') 
                    start_line = line_num + 1
                  end
                  if line_num > start_line
                    break if line == "\n"
                    begin
                      results = line.split
                      puts "Result0:"+results[0]
                      puts "Result1:"+results[1]
                      puts results[0].split('(')[1]
                      position_one= results[0].split('(')[1].gsub(')','').to_i - 1
                      aasequence_one = AAsequence.first(:seq_id=> seq.seq_id, :original_position=>position_one)
                      puts results[1].split('(')[1]
                      position_two = results[1].split('(')[1].gsub(')','').to_i - 1
                      aasequence_two = AAsequence.first(:seq_id=> seq.seq_id, :original_position=>position_two)
                      mean_one = results[2].to_f
                      mean_two = results[3].to_f
                      correlation = results[4].to_f
                      NewCap.create(:aasequence_id => aasequence_one.AAsequence_id,
                                  :position_one => position_one,
                                  :position_two => position_two,
                                  :mean_one => mean_one,
                                  :mean_two => mean_two,
                                  :correlation => correlation,
                                  :seq_id => seq.seq_id )
                    rescue Exception => e
                      puts e.message
                    end
                  end
                  line_num +=1
                end #end while
              end #end if
            end
        end
       }
    end
    thread_array.map{|t| t.join}
  end
  
  def calculate_disorder_consensus_threaded(thread_num=65, second_thread_num=10)
    alignment_array =[]
    Alignment.all(:alignment_name => self.alignment_name, :order => [:align_order.asc]).each do |alignment|
      alignment_array << alignment
    end
    thread_array=[]
    thread_num.times do |i|
        thread_array[i] = Thread.new{
          while alignment_array.length > 0 do
            alignment = alignment_array.pop
            puts alignment.sequence.abrev_name + ":STARTED"
            alignment.sequence.calculate_disorder_consensus_threaded(second_thread_num)
            puts alignment.sequence.abrev_name + ":DONE"
          end
        }
     end
     thread_array.map{|t| t.join}          
  end
  
  def calculate_intraresidue_consensus_threaded(thread_num=65,special=false)
    alignment_array =[]
    Alignment.all(:alignment_name => self.alignment_name, :order => [:align_order.asc]).each do |alignment|
      alignment_array << alignment
    end
    thread_array=[]
    thread_num.times do |i|
        thread_array[i] = Thread.new{
          while alignment_array.length > 0 do
            alignment = alignment_array.pop
            if PercentIdentity.all(:seq1_id => alignment.sequence.seq_id, :percent_id.gte => 19,:percent_id.lt => 90, :alignment_name => alignment.alignment_name).count > 9
              puts alignment.sequence.abrev_name + ":Started"
              alignment.sequence.calculate_intra_consensus(special)
              puts alignment.sequence.abrev_name + ":DONE"
            end
          end
        }
     end
     thread_array.map{|t| t.join}          
  end
  
  def evaluate_caps_distances(thread_num = 200)
    self.sequences.each do |s|
      s.evaluate_new_caps_distances(thread_num)
    end
  end

  def caps_report
    self.sequences.each do |s|
      caps_count = s.new_caps_count
      gtt = s.new_caps_count_gt_twenty
      gtf = s.new_caps_count_gt_50
      gth = s.new_caps_count_gt_100
      puts s.abrev_name + "\t| Caps Count: " + caps_count.to_s + "\t | Caps GTT Count: " + gtt.to_s + ", (%20): " + (caps_count > 0 ? (gtt/caps_count).to_s : "0") +  "\t|CAPS GTF: " + gtf.to_s + ", (%50): " + (caps_count > 0 ? (gtf/caps_count).to_s : "0") + "\t|CAPS GTH: " + gth.to_s + ", (%100): " + (caps_count > 0 ? (gth/caps_count).to_s : "0") 
    end
  end
  
  
  def import_conseq(thread_num=65)
    seq_array = []
    self.sequences.each do |s|
      seq_array << s
    end
    thread_array=[]
    thread_num.times do |i|
      thread_array[i] = Thread.new{
        while seq_array.length > 0 do
          seq = seq_array.pop      
          if PercentIdentity.all(:seq1_id => seq.seq_id, :percent_id.gte => 19,:percent_id.lt => 90, :alignment_name => self.alignment_name).count > 9
            #open file that corresponds to this sequence
            puts filename = "temp_data/#{self.alignment_name}/conseq/#{seq.abrev_name}.conseq"
            if File.exists?(filename)
              #seq.conseqs.destroy!
              Conseq.all(:seq_id=>seq.seq_id)
              puts "File exists"
              file = File.new(filename, "r")
              start_line = 13
              line_num = 1
              while (line = file.gets)
                if line_num > start_line
                  results = line.split     
                  puts line
                  begin       
                  cq = Conseq.new(:seq_id =>seq.seq_id,
                                :aasequence_id => seq.a_asequences.first(:original_position=>(results[0].to_i-1)).id,
                                :position=> results[0].to_i-1,
                                :score => results[2].to_f,
                                :color => results[3].to_i,
                                :state => results[4],
                                :function => results[5],
                                :msa_data => results[6],
                                :residue_variety => results[7])
                  cq.valid?
                  puts cq.errors.inspect()
                  cq.save
                  rescue Exception => e
                   puts seq.abrev_name + line
                   puts e.message
                   break  
                  end
                end
                line_num +=1
              end #end while
            end #end if
          end
        end
      }
    end
    thread_array.map{|t| t.join}
  end
  
  # def import_rate4site(thread_num=65)
  #   seq_array = []
  #   self.sequences.each do |s|
  #     seq_array << s
  #   end
  #   thread_array=[]
  #   thread_num.times do |i|
  #     thread_array[i] = Thread.new{
  #       while seq_array.length > 0 do
  #         seq = seq_array.pop      
  #         if PercentIdentity.all(:seq1_id => seq.seq_id, :percent_id.gte => 19,:percent_id.lt => 90, :alignment_name => self.alignment_name).count > 9
  #           #open file that corresponds to this sequence
  #           puts filename = "temp_data/#{self.alignment_name}/#{seq.abrev_name}.conseq"
  #           if File.exists?(filename)
  #             #seq.conseqs.destroy!
  #             Conseq.all(:seq_id=>seq.seq_id).destroy!
  #             puts "File exists"
  #             file = File.new(filename, "r")
  #             start_line = 12
  #             line_num = 1
  #             while (line = file.gets)
  #               if line_num > start_line
  #                 results = line.split     
  #                 puts line
  #                 begin       
  #                 cq = Conseq.new(:seq_id =>seq.seq_id,
  #                               :aasequence_id => seq.a_asequences.first(:original_position=>(results[0].to_i-1)).id,
  #                               :position=> results[0].to_i-1,
  #                               :score => results[2].to_f,
  #                               :msa_data => results[5],
  #                               :std =>results[4].to_f,
  #                               :qq_interval =>results[3])
  #                 cq.valid?
  #                 puts cq.errors.inspect()
  #                 cq.save
  #                 rescue Exception => e
  #                  puts seq.abrev_name + line
  #                  puts e.message
  #                  break  
  #                 end
  #               end
  #               line_num +=1
  #             end #end while
  #           end #end if
  #         end
  #       end
  #     }
  #   end
  #   thread_array.map{|t| t.join}
  # end
  
  def import_cornet(thread_num=65)
    seq_array = self.sequences.to_a
    thread_array=[]
    thread_num.times do |i|
      thread_array[i] = Thread.new{
        while seq_array.length > 0 do
          seq = seq_array.pop      
          if PercentIdentity.all(:seq1_id => seq.seq_id, :percent_id.gte => 19,:percent_id.lt => 90, :alignment_name => self.alignment_name).count > 9
            #open file that corresponds to this sequence
            puts filename = "temp_data/#{self.alignment_name}/cornet/#{seq.abrev_name}.cornet"
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
        end
      }
    end
    thread_array.map{|t| t.join}
  end
  
  def generate_jalview_annotation_consensus()
    File.open("#{self.alignment_name}.gff", 'w') do |f1| 
      f1.puts "both\t00FF00"
      f1.puts "cicp\t66FFCC" 
      f1.puts "low_disorder\tFFFF00" 
      f1.puts "avg_disorder\tFFCC00" 
      f1.puts "medium_disorder\tFF9900"
      f1.puts "highly_disordered\tFF6600" 
      f1.puts "extremely_disordered\tFF0000"
      self.sequences.each do |seq|     
        seq.a_asequences.all.each do |aa|
          # if aa.disorder_consensus > 0.5 && aa.disorder_consensus < 0.6
          #   feature_type = "low_disorder"
          # elsif aa.disorder_consensus >= 0.6 && aa.disorder_consensus < 0.7
          #   feature_type = "avg_disorder"
          # elsif aa.disorder_consensus >= 0.7 && aa.disorder_consensus < 0.8
          #   feature_type = "medium_disorder"
          # elsif aa.disorder_consensus >= 0.8 && aa.disorder_consensus < 0.9
          #   feature_type = "highly_disordered"
          # elsif aa.disorder_consensus > 0.9
          #   feature_type = "extremely_disordered"
          # else
          #   feature_type = "no_disorder"
          # end
          feature = aa.jalview_residue_color
          if feature != "no_disorder"
            f1.puts "Jalview\t#{seq.abrev_name}\t-1\t#{aa.original_position}\t#{aa.original_position}\t#{feature}"
          end
       end
     end
   end
   #return jalview_string
  end


  def fix_alignment
    alignments = Alignment.all(:alignment_name => self.alignment_name)
    f = File.new(self.alignment_name+".fa", "w")
    @pos_file = File.new(self.alignment_name+"_added_positions.txt", "w")
    alignments.each do |a|
      @pos_file.write("****#{a.sequence.abrev_name}\n")
      seq = a.sequence.sequence.gsub("\n",'').gsub("\r",'')
      aaseq = a.alignment_sequence
      new_alignment = ""
      seq_pos = 0
      seq.each_char do |s|
        while aaseq[seq_pos] == "-"
          new_alignment << aaseq[seq_pos]
          seq_pos +=1
        end
        if !aaseq[seq_pos].nil?
          if s.downcase == aaseq[seq_pos].downcase
            new_alignment << aaseq[seq_pos]
            seq_pos +=1
          else
            new_alignment << s
            @pos_file.write(s +":#{seq_pos}\n")
          end
        end
      end
      f.write(">"+a.sequence.abrev_name+"\n")
      f.write(new_alignment+"\n")
    end
    f.close
    @pos_file.close
  end
  
  @queue= :disicc
  #resque task
  def self.perform(alignment_id)
    alignment = Alignment.get(alignment_id)
    alignment.run_and_store_disorder(thread_num=20)
    alignment.run_align_assess
    alignment.run_xdet_threaded(dir="temp_data", thread_num=4)
    alignment.import_xdet_new(thread_num=4)
    alignment.run_rate4site
    alignment.import_rate4site(thread_num=4)
    alignment.run_perl_caps
    alignment.import_caps_threaded(dir = "temp_data/#{self.alignment_name}/", thread_num=4)
  end

end
