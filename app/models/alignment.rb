class Alignment
  include DataMapper::Resource
  
  property :align_id, Serial
  property :seq_id, Integer, :required => true
  property :alignment_name, String, :required => true
  property :align_order, Integer, :required => true
  property :alignment_sequence, Text, :required => true, :default => ""
  property :fasta_title, Text, :required => false, :default => ""
  property :deleted_at, ParanoidDateTime
  
  # has n, :disorder_values
  has 1, :sequence, 'Sequence', :child_key => [:seq_id]
  
  #validates_uniqueness_of :alignment_name
  
  #runs the score app to get percent identities for each sequence
  def run_score
    filename = self.generate_fasta_alignment_file
    string = "./lib/score_mac #{filename} temp_data/#{self.alignment_name}_res.txt temp_data/#{self.alignment_name}_dif.txt temp_data/#{self.alignment_name}_alignments.txt"
    puts string
    if system(string)
      
    end
  end
  
  def sequence
    Sequence.get(self.seq_id)
  end
  
  def sequences
    alignments = Alignment.all(:alignment_name => self.alignment_name, :order=>[:align_order])
    seq_ids = alignments.map{|a| a.seq_id}
    Sequence.all(:id=>seq_ids)
  end
  
  
  #runs the disorder apps for every sequence in the alignment
  #then calculates the disorder consequence
  # @params id
  def run_consensus_disorder
    self.sequences.each do |sequence|
        sequence.run_and_store_disorder()
        sequence.calculate_disorder_consensus()
    end
  end
  
  #runs the AlignAssess app to get percent identities for each sequence in the alignment
  def run_align_assess
    filename = self.generate_fasta_alignment_file_for_all
    string = "./lib/AlignAssess_wShorterID #{filename} P"
    seq_array = Array.new
    if system(string)
      seq_id_array = self.sequences.map{|s| s.id}
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
        for i in 1..range#(row_num)          
          PercentIdentity.first_or_create(:seq1_id=>seq_id_array[row_num],
                                                    :seq2_id=>seq_id_array[i],
                                                    :alignment_name => self.alignment_name,
                                                    :percent_id=>seq_array[row_num][i])
          # print "[#{row_num}:#{i-1}=>#{row[i]}],"
        end
        #print "\n"
      end
    end
  end
  def run_caps_mac
    self.run_align_assess
    Dir.mkdir("temp_data/#{self.alignment_name}") unless File.directory?("temp_data/#{self.alignment_name}")
    alignments = Alignment.all(:alignment_name => self.alignment_name)
    alignments.each do |alignment|
      alignment.generate_pid_fasta_file("temp_data/#{self.alignment_name}")
    end
    system "./lib/comp_apps/caps/caps_mac -F temp_data/#{self.alignment_name} --intra"
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
            system "./lib/comp_apps/caps2/caps2 -F temp_data/#{self.alignment_name}_#{this_seq_id} --intra"
            puts "Finsihed Caps2 for #{this_seq_id}"
          end
       }
     end
     thread_array.map{|t| t.join}
    end
  
  def run_xdet
    self.run_align_assess
    Dir.mkdir("temp_data/#{self.alignment_name}") unless File.directory?("temp_data/#{self.alignment_name}")
    alignments = Alignment.all(:alignment_name => self.alignment_name)
    alignments.each do |alignment|
      filename= alignment.generate_pid_fasta_file("temp_data/#{self.alignment_name}")
      system "./lib/comp_apps/XDet/xdet_linux32 #{filename} ~/Rails/DisICC/lib/comp_apps/XDet/Maxhom_McLachlan.metric >> #{filename}_xdet"
    end
  end
  
  def run_rate4site
    self.run_align_assess
    Dir.mkdir("temp_data/#{self.alignment_name}") unless File.directory?("temp_data/#{self.alignment_name}")
    alignments = Alignment.all(:alignment_name => self.alignment_name)
    alignments.each do |alignment|
      filename= alignment.generate_pid_fasta_file("temp_data/#{self.alignment_name}")
      system "./lib/comp_apps/conseq/rate4site_fast -s #{filename} -o #{filename}_conseq"
    end
  end
  
  def generate_pid_fasta_files
    self.run_align_assess
    self.sequences.each do |seq|
      fasta_string=""
      pids = PercentIdentity.all(:seq1_id => seq.seq_id, :percent_id.gte => 20, :order =>[:percent_id.desc])
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
  
  #assumes that run_align_asses has already occured
  def generate_pid_fasta_file(dir="temp_data", longest_alignment_length=0)
    fasta_string=""
    seq = Sequence.get(self.seq_id)
    pids = PercentIdentity.all(:seq1_id => self.seq_id, :percent_id.gte => 20, :order =>[:percent_id.desc],:unique=>true)
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
  
  def generate_fasta_alignment_file_for_all(filename="")
    alignments = Alignment.all(:alignment_name=> self.alignment_name)
    if filename.empty?
      filepath = "temp_data/"+self.alignment_name+"_alignment"+Time.now.to_i.to_s+".fasta"
    else
      filepath = "temp_data/#{filename}"
    end
    f = File.new(filepath, "w+")
    alignments.each do |alignment|
      f.write(alignment.fasta_alignment_string)
    end
    f.close
    filepath
  end
  

  def import_caps
    #find correct directory
    #dir_name = "#{root_dir}/#{self.alignment_name}"
    #for each sequence in the alignment
    self.sequences.each do |seq|
      #open file that corresponds to this sequence
      puts filename = "#{self.alignment_name}_#{seq.abrev_name}_pid.fasta.out"
      if File.exists?(filename)
        puts "File exists"
        file = File.new(filename, "r")
        start_line = 99999999999
        line_num = 1
        while (line = file.gets)
          if line.include?('Position')
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
            NewCap.create(:aasequence_id => aasequence_one.AAsequence_id,
                        :position_one => position_two,
                        :position_two => position_one,
                        :mean_one => mean_one,
                        :mean_two => mean_two,
                        :correlation => correlation,
                        :seq_id => seq.seq_id )
          end
          line_num +=1
        end #end while
      end #end if
    end #end sequences.each
  end


end
