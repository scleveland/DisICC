class AlignmentsController < ApplicationController
  # GET /alignments
  # GET /alignments.xml
  def index
    @alignments = Alignment.all#(:seq_id=>current_user.sequences.map{|s| s.seq_id})
    sql = "Select Distinct alignment1_id, alignment2_id From inter_caps WHERE alignment1_id IS NOT NULL"
    results = repository.adapter.select(sql)
    @al1_hash = {}
    results.each do |r|
      puts "#{r.alignment1_id}, #{r.alignment2_id}"
      @al1_hash["#{Alignment.get(r.alignment1_id.to_i).alignment_name}"] =[]
      @al1_hash["#{Alignment.get(r.alignment2_id.to_i).alignment_name}"] =[]
    end
    results.each do |r|
      if !r.alignment1_id.nil? && !r.alignment2_id.nil?
        @al1_hash["#{Alignment.get(r.alignment1_id.to_i).alignment_name}"] << Alignment.get(r.alignment2_id.to_i).alignment_name
        @al1_hash["#{Alignment.get(r.alignment2_id.to_i).alignment_name}"] << Alignment.get(r.alignment1_id.to_i).alignment_name
        
      end
    end
    @al1_hash.each do |al|
      al.uniq!
    end
    respond_to do |format|
      format.html # index.html.erb
      format.xml  { render :xml => @alignments }
    end
  end

  # GET /alignments/1
  # GET /alignments/1.xml
  def show
    @alignment = Alignment.find(params[:id])

    respond_to do |format|
      format.html # show.html.erb
      format.xml  { render :xml => @alignment }
    end
  end

  # GET /alignments/new
  # GET /alignments/new.xml
  def new
    @alignment = Alignment.new
    @seq_types = {:N=>"N", :L=>"L", :P => "P", "P Cut"=> "P Cut"}
    respond_to do |format|
      format.html # new.html.erb
      format.xml  { render :xml => @alignment }
    end
  end

  # GET /alignments/1/edit
  def edit
    @alignment = Alignment.find(params[:id])
  end

  # POST /alignments
  # POST /alignments.xml
  def create
    @alignment = Alignment.new(params[:alignment])

    respond_to do |format|
      if @alignment.save
        format.html { redirect_to(@alignment, :notice => 'Alignment was successfully created.') }
        format.xml  { render :xml => @alignment, :status => :created, :location => @alignment }
      else
        format.html { render :action => "new" }
        format.xml  { render :xml => @alignment.errors, :status => :unprocessable_entity }
      end
    end
  end

  # PUT /alignments/1
  # PUT /alignments/1.xml
  def update
    @alignment = Alignment.find(params[:id])

    respond_to do |format|
      if @alignment.update_attributes(params[:alignment])
        format.html { redirect_to(@alignment, :notice => 'Alignment was successfully updated.') }
        format.xml  { head :ok }
      else
        format.html { render :action => "edit" }
        format.xml  { render :xml => @alignment.errors, :status => :unprocessable_entity }
      end
    end
  end

  # DELETE /alignments/1
  # DELETE /alignments/1.xml
  def destroy
    @alignment = Alignment.get(params[:id])
    alignments = Alignment.all(:alignment_name => @alignment.alignment_name)
    alignments.each do |a|
      a.deleted_at = Time.now
      a.save
    end

    respond_to do |format|
      format.html { redirect_to(alignments_url) }
      format.xml  { head :ok }
    end
  end
  
  def run_disorder
    Alignment.get(params[:id]).run_consensus_disorder
    redirect_to(display_disorder_annotated_alignment_alignment_path(params[:id]))
  end
  
  def residue_color(dis_avg, con_avg)
    if con_avg > 0.75
      if dis_avg >= 0.5
       @color = "00FF00"  #3394560  3CF3C0
      else
       @color =  "66FFCC" #6750156  66FFCC
      end
    else  #color for disorder only
      if dis_avg >= 0.5 && dis_avg < 0.6
       @color = "FFFF00" #16776960  FFFF00
      elsif dis_avg >= 0.6 && dis_avg < 0.7
       @color =  "FFCC00" #16763904  FFCC00
      elsif dis_avg >= 0.7 && dis_avg < 0.8
       @color =  "FF9900" #16750848  FF9900
      elsif dis_avg >= 0.8 && dis_avg < 0.9
       @color = "FF6600"  #16737792  FF6600
      elsif dis_avg >= 0.9
       @color =  "FF0000" #16711680  FF0000
      else
       @color = "999999"
      end
    end
    return @color
  end
  
  def alignment_to_positions(alignment)
    counter = 0
    alignment.alignment_sequence.each_char do |aa|
      if aa != "-"
        aaseq = AAsequence.first(:seq_id=>alignment.seq_id, :original_position=> counter)
        if !aaseq.nil?
          AlignmentPosition.create(:alignment_id => alignment.align_id,
                            :position => counter,
                            :aasequence_id => aaseq.AAsequence_id)
        end
      end
      counter +=1
    end
  end
  
  def upload
    @seq_types = {:N=>"N", :L=>"L", :P => "P"}
  end
  
  def disorder_brief_report
    thread_num=70
    @dis_array=[]
    alignment_array=[]
    Alignment.all(:alignment_name => Alignment.get(params[:id]).alignment_name, :order => [:align_order.asc]).each do |alignment|
      alignment_array << alignment
    end
    thread_array=[]
    thread_num.times do |i|
        thread_array[i] = Thread.new{
          while alignment_array.length > 0 do
            alignment = alignment_array.pop  
            temp_hash={}
            seq = alignment.sequence
            seq.disorders.each do |disorder|
              temp_hash[disorder.disorder_type] = repository(:default).adapter.select("SELECT COUNT (*) FROM disorder_values where disorder_id =#{disorder.id}")#disorder.disorder_values.count
              temp_hash["#{disorder.disorder_type}_id"] =disorder.disorder_id
            end
            temp_hash[:name] =seq.abrev_name
            temp_hash[:id] = seq.id
            temp_hash[:type] = seq.seq_type
            temp_hash[:align_id] = alignment.align_id
            temp_hash[:consensus_count] = seq.a_asequences.all(:disorder_consensus.gte => 0.5).count
            temp_hash[:disorder_count] = seq.disorders.count
            @dis_array[alignment.align_order] = temp_hash
         end
        }
    end
    thread_array.map{|t| t.join}
    @disorder_types = repository(:default).adapter.select('SELECT DISTINCT disorder_type FROM Disorders')
    #Disorder.all(:fields=>[:disorder_type, :id], :unique=>true).map{|d| d.disorder_type}
  end
  
  def compensatory_brief_report
    thread_num=70
    @comp_array=[]
    alignment_array=[]
    @alignment_name = Alignment.get(params[:id]).alignment_name
    Alignment.all(:alignment_name => @alignment_name, :order => [:align_order.asc]).each do |alignment|
      alignment_array << alignment
    end
    thread_array=[]
    thread_num.times do |i|
        thread_array[i] = Thread.new{
          while alignment_array.length > 0 do
            alignment = alignment_array.pop  
            temp_hash={}
            seq = alignment.sequence
            puts seq.abrev_name + ":intra"
            temp_hash[:intra] = repository(:default).adapter.select("SELECT COUNT (*) FROM intra_residue_contacts where seq_id =#{seq.id}")#seq.intra_residue_contacts.count
            puts seq.abrev_name + ":new_caps"
            temp_hash[:new_caps] = repository(:default).adapter.select("SELECT COUNT (*) FROM new_caps where seq_id =#{seq.id}")#seq.new_caps.count
            puts seq.abrev_name + ":xdet"
            temp_hash[:xdet] = repository(:default).adapter.select("SELECT COUNT (*) FROM xdets where seq_id =#{seq.id}")#seq.xdets.count
            puts seq.abrev_name + ":conseq"
            temp_hash[:conseq] = repository(:default).adapter.select("SELECT COUNT (*) FROM conseqs where seq_id =#{seq.id}")#seq.conseqs.count
            temp_hash[:name] = seq.abrev_name
            temp_hash[:id] = seq.seq_id
            temp_hash[:pids] = PercentIdentity.all(:seq1_id => alignment.sequence.seq_id, :percent_id.gt => 19,:percent_id.lt => 90, :alignment_name => alignment.alignment_name).count
            temp_hash[:aa_count] = seq.a_asequences.count
            temp_hash[:seq_length] = seq.sequence.gsub("\n","").gsub("\r","").length
            temp_hash[:alignment_position_count] = AlignmentPosition.all(:alignment_id => alignment.align_id).count
            temp_hash[:consensus] = alignment.sequence.a_asequences.all(:contact_consensus.gte => 0.5).count
            @comp_array[alignment.align_order] = temp_hash
         end
        }
    end
    thread_array.map{|t| t.join}
  end
  
  def calculate_pids
    alignment = Alignment.get(params[:id])
    alignment.run_align_assess
    redirect_to(percent_identities_alignment_path(params[:id]))
  end
  
  def pre_process_fasta_file
    @name = Time.now.to_s + params[:datafile].original_filename
    directory = "temp_data"
    @new_file = File.join(directory,@name)
    @alignment_name = params[:alignment_name]
    @seq_type = params[:seq_type]
  
      sequences = Sequence.all(:seq_type => @seq_type)
    
    @seq_options = Hash.new
    sequences.map{|k| @seq_options = @seq_options.merge({k.seq_name => k.id.to_s})}
    File.open(@new_file, "wb"){ |f| f.write(params['datafile'].read)}
    logger.debug "HELLYEAH"
    begin
      file = File.new(@new_file, "r")
      @fasta_name_arrays = Array.new
      while (line = file.gets)
        if line.count(">") > 0
          @fasta_name_arrays << line.gsub(">", "")
        end
      end
      file.close
      rescue => err
          puts "Exception: #{err}"
          err
    end
  end
  
  
  def complete_process_fasta_file
    directory = "temp_data"
    new_file = File.join(directory,params[:datafile_name])
    fasta_hash = Hash.new
    range = params[:seq_num].to_i - 1
    (0..range).each do |i|
      fasta_hash = fasta_hash.merge(params["fasta_name"+i.to_s].strip => params["seq"+i.to_s])
    end
    # logger.debug "BEFORE*************************"
    # begin
    #   file = File.new(@new_file, "r")
    #   logger.debug "After FILE OPEN"
    #   puts "After FILE OPEN"
    #   counter = 0
    #   order_count = 0
    #   abrev_name = ""
    #   alignment_sequence = ""
    #   fasta_hash = Hash.new
    #   logger.debug "BEFORE FASTA HASH*****"
    #   puts "BEFORE FASTA HASH*****"
    #   logger.debug "BEFORE@1232132131****"
    #   range = params[:seq_num].to_i - 1
    #   logger.debug "BEFORE FASTA RANGE***********"
    #   puts "BEFORE FASTA RANGE***********"
    #   (0..range).each do |i|
    #     logger.debug "BLAH"
    #     fasta_hash = fasta_hash.merge(params["fasta_name"+i.to_s].strip => params["seq"+i.to_s])
    #   end
    #   logger.debug {fasta_hash}
    #   logger.debug "AFTER FASTA HASH*******"
    #   while (line = file.gets)
    #     if line.count(">") > 0 && counter > 0
    #       #save the current sequence to an alignment
    #       logger.debug {fasta_hash[abrev_name]}
    #       puts fasta_hash[abrev_name]
    #       @sequence = Sequence.get(fasta_hash[abrev_name.strip])
    #       logger.debug "After"
    #       puts "After"  
    #       logger.debug "OHNO NONONONONO" 
    #       @alignment = Alignment.new(:seq_id => @sequence.id,
    #                        :alignment_name => params[:alignment_name],
    #                        :align_order => order_count,
    #                        :alignment_sequence => alignment_sequence,
    #                        :fasta_title => abrev_name)
    #        logger.debug "SHAZZZZZAAAAAAAM"                 
    #       logger.debug { @alignment.valid?}
    #       if !@alignment.valid?
    #         puts @alignment.errors.inspect()
    #       end
    #       logger.debug "VALID"
    #       logger.debug { @alignment.errors.inspect }
    #       @alignment.save
    #       alignment_to_positions(@alignment)              
    #       #this is the sequene label
    #       abrev_name = line.gsub(">", "")
    # 
    #       order_count += 1
    #       alignment_sequence = ""
    #     elsif line.count(">") > 0
    #       abrev_name = line.gsub(">","")
    #       logger.debug { abrev_name }
    #       alignment_sequence =""
    #     elsif counter > 0
    #       alignment_sequence = alignment_sequence + line.lstrip.rstrip
    #     end
    #     counter = counter + 1
    #   end
    #    puts fasta_hash[abrev_name]
    #     @sequence = Sequence.get(fasta_hash[abrev_name.strip])
    #     logger.debug "After"
    #     puts "After"  
    #     logger.debug "OHNO NONONONONO" 
    #     @alignment = Alignment.new(:seq_id => @sequence.seq_id,
    #                      :alignment_name => params[:alignment_name],
    #                      :align_order => order_count,
    #                      :alignment_sequence => alignment_sequence,
    #                      :fasta_title => abrev_name)
    #      logger.debug "SHAZZZZZAAAAAAAM"                 
    #     logger.debug { @alignment.valid?}
    #     if !@alignment.valid?
    #       puts @alignment.errors.inspect()
    #     end
    #     logger.debug "VALID"
    #     logger.debug { @alignment.errors.inspect }
    #     @alignment.save
    #     @alignment.alignment_to_positions              
    #     #this is the sequene label
    #     abrev_name = line.gsub(">", "")
    # 
    #     order_count += 1
    #     alignment_sequence = ""
    #   file.close
    #   rescue => err
    #       puts "Exception: #{err}"
    #       err
    # end
    require 'process_fasta.rb'
    Resque.enqueue(ProcessFasta, new_file, fasta_hash, params[:alignment_name], current_user.id)
    flash[:notice] = "You will recieve an email when you Alignment is ready."
    redirect_to(alignments_path)
  end
  
  def process_fasta_file_and_save_sequences
    name = Time.now.to_s + params[:datafile].original_filename
    directory = "temp_data"
    new_file = File.join(directory,name)
    File.open(new_file, "wb"){ |f| f.write(params['datafile'].read)}
    #sequence_list = Sequence.create_sequences_from_fasta(new_file, current_user.id)
    #sequences = Sequence.all(:seq_id => sequence_list.map{|s| s.seq_id})
    #sequences = []
    file = File.new(new_file, 'r')
    ff = Bio::FlatFile.new(Bio::FastaFormat, file)
    order_count = 0
    ff.each_entry do |f|
      puts f.definition
      seq = Sequence.create_sequence_from_bioruby_fasta_entry(f, params[:seq_type],current_user.id)
      puts seq.abrev_name
      alignment = Alignment.create(:seq_id => seq.id,
                        :alignment_name => params[:alignment_name],
                        :align_order => order_count,
                        :alignment_sequence =>f.naseq)
      #alignment_to_positions(alignment) 
      alignment.alignment_to_positions  
      order_count += 1
      begin
        seq.run_and_store_disorder()
        
      rescue Exception => e  
        puts self.abrev_name + " Diorder Generation Failed"
        puts e.message
      end
      begin
        seq.calculate_disorder_consensus()
      rescue Exception => e  
        puts self.abrev_name + " Diorder Consensus Calculation Failed"
        puts e.message
      end
    end
    redirect_to(alignments_path)
  end
  
  def process_fasta_file
    name = Time.now.to_s + params[:datafile_name]#.original_filename
    directory = "temp_data"
    @new_file = File.join(directory,params[:datafile_name])
    File.open(@new_file, "wb"){ |f| f.write(params['datafile_name'])}
    logger.debug "HELLYEAH"
    begin
      file = File.new(@new_file, "r")
      counter = 0
      order_count = 0
      abrev_name = ""
      alignment_sequence = ""
      while (line = file.gets)
        if line.count(">") > 0 && counter > 0
          #save the current sequence to an alignment
          logger.debug { alignment_sequence }
          logger.debug "Before"
          if abrev_name.count("/") > 0
            abrev_name = abrev_name.split("/")[0]
          end
          @sequence = Sequence.first(:abrev_name => abrev_name.lstrip.rstrip, :seq_type => params[:seq_type])
          logger.debug "After"
          logger.debug { @sequence.to_s}
          logger.debug "OHNO" 
          @alignment = Alignment.new(:seq_id => @sequence.id,
                           :alignment_name => params[:alignment_name],
                           :align_order => order_count,
                           :alignment_sequence => alignment_sequence)
          @alignment.valid?
          logger.debug "VALID"
          logger.debug { @alignment.errors.inspect }
          @alignment.save
          #alignment_to_positions(@alignment)
          @alignment.alignment_to_positions            
          #this is the sequene label
          abrev_name = line.gsub(">", "")

          order_count += 1
          alignment_sequence = ""
        elsif line.count(">") > 0
          abrev_name = line.gsub(">","")
          logger.debug { abrev_name }
          alignment_sequence =""
        elsif counter > 0
          alignment_sequence = alignment_sequence + line.lstrip.rstrip
        end
        counter = counter + 1
      end
      file.close
      rescue => err
          puts "Exception: #{err}"
          err
    end
    redirect_to(alignments_path)
  end
  
  
  def calculate_intraresidue_consensus
    Alignment.all(:alignment_name => Alignment.get(params[:id]).alignment_name, :order => [:align_order.asc]).each do |alignment|
      AAsequence.all(:seq_id => alignment.seq_id).each do |aaseq|
        count = 0
        if !IntraResidueContact.first(:seq_id => aaseq.seq_id, :first_residue=> aaseq.original_position).nil?
          count +=1
        elsif !IntraResidueContact.first(:seq_id => aaseq.seq_id, :second_residue=> aaseq.original_position).nil?
        end
        if !Conseq.first(:aasequence_id => aaseq.AAsequence_id).nil?
          if Conseq.first(:aasequence_id => aaseq.AAsequence_id).color < 4
            count +=1
          end
        end
        if !Xdet.first(:aasequence_id => aaseq.AAsequence_id).nil?
          if Xdet.first(:aasequence_id => aaseq.AAsequence_id).correlation > 0.0 || Xdet.first(:aasequence_id => aaseq.AAsequence_id).correlation == -2
            count +=1
          end
        end
        if !NewCap.first(:seq_id=> aaseq.seq_id, :position_one => aaseq.original_position).nil?
          count +=1
        elsif !NewCap.first(:seq_id=> aaseq.seq_id, :position_two => aaseq.original_position).nil?
          count +=1
        end
        aaseq.contact_positive_consensus = count /4
        aaseq.save
      end  
    end  
    redirect_to(alignments_path)                        
  end  
  
  
  def calculate_intraresidue_consensus_threaded(thread_num=65)
    alignment_array =[]
    Alignment.all(:alignment_name => Alignment.get(params[:id]).alignment_name, :order => [:align_order.asc]).each do |alignment|
      alignment_array << alignment
    end
    thread_array=[]
    thread_num.times do |i|
        thread_array[i] = Thread.new{
          while alignment_array.length > 0 do
            alignment = alignment_array.pop

            alignment.sequence.calculate_intra_consensus
            # puts alignment.sequence.abrev_name + ":STARTED"
            #             AAsequence.all(:seq_id => alignment.seq_id).each do |aaseq|
            #               count = 0
            #               # if !IntraResidueContact.first(:seq_id => aaseq.seq_id, :first_residue=> aaseq.original_position).nil?
            #               #                 count +=1
            #               #               elsif !IntraResidueContact.first(:seq_id => aaseq.seq_id, :second_residue=> aaseq.original_position).nil?
            #               #                 count +=1
            #               #               end
            #               if !Conseq.first(:aasequence_id => aaseq.AAsequence_id).nil? && Conseq.first(:aasequence_id => aaseq.AAsequence_id).color > 4
            #                   count +=1
            #               end
            #               if !Xdet.first(:aasequence_id => aaseq.AAsequence_id).nil? & (Xdet.first(:aasequence_id => aaseq.AAsequence_id).correlation > 0.0 || Xdet.first(:aasequence_id => aaseq.AAsequence_id).correlation == -2)
            #                   count +=1
            #               end
            #               if !NewCap.first(:seq_id=> aaseq.seq_id, :position_one => aaseq.original_position).nil? || !NewCap.first(:seq_id=> aaseq.seq_id, :position_two => aaseq.original_position).nil?
            #                 count +=1
            #               end
            #               aaseq.contact_consensus = count /3#4
            #               aaseq.save
            #end  

            puts alignment.sequence.abrev_name + ":DONE"
          end
        }
     end
     thread_array.map{|t| t.join}          
     redirect_to(alignments_path)   
  end
  
  def percent_identities 
    @align = Alignment.get(params[:id])
    @pids = {}
    @alignments =  Alignment.all(:alignment_name => @align.alignment_name, :order=>[:align_order])
    @alignments.each do |alignment|
      @pids[alignment.sequence.id] =  PercentIdentity.all(:seq1_id => alignment.sequence.id, :percent_id.gt => 19,:percent_id.lt => 100, :alignment_name => @align.alignment_name)
    end
  end
  
  def calculate_disorder_consensus_threaded
    thread_num = 65
    alignment_array =[]
    Alignment.all(:alignment_name => Alignment.get(params[:id]).alignment_name, :order => [:align_order.asc]).each do |alignment|
      alignment_array << alignment
    end
    thread_array=[]
    thread_num.times do |i|
        thread_array[i] = Thread.new{
          while alignment_array.length > 0 do
            alignment = alignment_array.pop
            puts alignment.sequence.abrev_name + ":STARTED"
            alignment.sequence.calculate_disorder_consensus_threaded(20)
            puts alignment.sequence.abrev_name + ":DONE"
          end
        }
     end
     thread_array.map{|t| t.join}          
     redirect_to(alignments_path)   
  end
   # Alignment.all(:alignment_name => Alignment.get(826).alignment_name, 
   #                              :order => [:align_order.asc]).each do |alignment|
   #    
   #    if !AAsequence.first(:seq_id => alignment.seq_id, :contact_consensus.gte => 0.75).nil?
   #      @seq_contact_count += 1
   #    end
   #    puts Sequence.first(:seq_id => alignment.seq_id).abrev_name + ":" + @seq_contact_count
   #  end
  
   def display_disorder_annotated_alignment
     thread_num = 65
     @display_array = Array.new
     @max_count = 0
     @contact_consensus_array = Array.new
     @seq_contact_count = Alignment.all(:alignment_name => Alignment.get(params[:id]).alignment_name).count
     longest_alignment = 0;
     alignment_array = []
     @alignment_name = Alignment.get(params[:id]).alignment_name
     Alignment.all(:alignment_name => Alignment.get(params[:id]).alignment_name, 
                                 :order => [:align_order.asc]).each do |alignment|
      puts alignment.alignment_sequence.length
      if alignment.alignment_sequence.length > longest_alignment
        longest_alignment = alignment.alignment_sequence.length
      end
      alignment_array << alignment
    end
    for i in 0..longest_alignment+1
      @contact_consensus_array[i] = Array.new(@seq_contact_count, 0)
    end
    #@contact_consensus_array = Array.new(longest_alignment, Array.new(@seq_contact_count,0))
    puts @contact_consensus_array.length
    puts "Into The Threads"
     thread_array=[]
      thread_num.times do |i|
        thread_array[i] = Thread.new{
          while alignment_array.length > 0 do
            alignment = alignment_array.pop
            sequence= alignment.sequence
            display_hash = Hash.new
            alignment_color_array = Array.new      
            cur_position = 0   
            orig_position = 0
            AlignmentPosition.all(:alignment_id => alignment.align_id, 
                         :order => [:alignment_position_id.asc]).each do |position|
             if position.position == cur_position
                amino_acid = sequence.a_asequences.first(:original_position=>orig_position) #AAsequence.first(:id => position.aasequence_id)
                unless amino_acid.nil?
                  alignment_color_array[cur_position] = residue_color(amino_acid.disorder_consensus, 0)
                  if @contact_consensus_array[cur_position][alignment.align_order].nil?
                    @contact_consensus_array[cur_position][alignment.align_order] = 0
                  end
                  if amino_acid.disorder_consensus >= 0.5
                    @contact_consensus_array[cur_position][alignment.align_order] = @contact_consensus_array[cur_position][alignment.align_order] + 1
                  end
                else
                  puts "Amino Acid doesn't exits: #{sequence.abrev_name} | #{cur_position}:#{orig_position}" 
                  alignment_color_array[cur_position] = residue_color(0, 0)
                  @contact_consensus_array[cur_position][alignment.align_order] = 0
                end
             else
                while position.position > cur_position
                             alignment_color_array[cur_position] = "FFFFFF"
                             cur_position += 1
                end
                amino_acid = sequence.a_asequences.first(:original_position=>orig_position) #AAsequence.first(:id => position.aasequence_id)
                unless amino_acid.nil?
                  alignment_color_array[cur_position] = residue_color(amino_acid.disorder_consensus, 0)
                  if @contact_consensus_array[cur_position].nil?
                    puts "OH no " + alignment.sequence.abrev_name
                  end
                  if @contact_consensus_array[cur_position][alignment.align_order].nil?
                     @contact_consensus_array[cur_position][alignment.align_order] = 0
                  end
                  if amino_acid.disorder_consensus >= 0.5
                     @contact_consensus_array[cur_position][alignment.align_order] = @contact_consensus_array[cur_position][alignment.align_order] + 1
                  end
                else
                    puts "Amino Acid doesn't exits: #{sequence.abrev_name} | #{cur_position}:#{orig_position}" 
                    alignment_color_array[cur_position] = residue_color(0, 0)
                    @contact_consensus_array[cur_position][alignment.align_order] = 0
                end
              end
              cur_position += 1
              orig_position +=1
              end                 
              puts display_hash["name"] = Sequence.first(:seq_id => alignment.seq_id).abrev_name 
              display_hash["alignment"] = alignment_color_array
              @display_array[alignment.align_order] = display_hash
            if @max_count < cur_position
                   @max_count = cur_position
            end
          end
        }
     end
     thread_array.map{|t| t.join}

     @contact_consensus_array = @contact_consensus_array.map{|a| a.inject(0){|sum,item| sum + item}}
     @cur_position = 0
     @tick_counter = 0
     @alignment_tick_array = Array.new
     while @cur_position <= @max_count
       @cur_position += 1
       @tick_counter += 1
       if @tick_counter != 25
         @alignment_tick_array << "FFFFFF"
       else
         @alignment_tick_array << "000000"
         @tick_counter = 0
       end
     end
     @display_hash = Hash.new
     @display_hash["name"] = ""
     @display_hash["alignment"] = @alignment_tick_array  
     @display_array << @display_hash
     if params[:aa_length].nil?
       @aa_length = 400
     else
       @aa_length = params[:aa_length].to_i
     end
     @ranges = (@max_count/@aa_length)

   end 

   def display_disorder_and_cicp_annotated_alignment
     thread_num = 65
     @display_array = Array.new
     @max_count = 0
     @cicp_and_disorder_stats = []
     @contact_consensus_array = Array.new
     @cicp_array = Array.new
     @cicp_contact_count =0
     align = Alignment.get(params[:id])
     @seq_contact_count = Alignment.all(:alignment_name =>align.alignment_name).count
     longest_alignment = 0;
     alignment_array = []
     unless Conservation.first(:alignment_name => align.alignment_name)
       align.run_aacon
     end
     if Conservation.first(:alignment_name => align.alignment_name)
       @conservation = Conservation.first(:alignment_name => align.alignment_name).results_array
     else
       @conservation = []
     end
     @alignment_name = Alignment.get(params[:id]).alignment_name
     Alignment.all(:alignment_name => Alignment.get(params[:id]).alignment_name, 
                                 :order => [:align_order.asc]).each do |alignment|
      puts alignment.alignment_sequence.length
      if alignment.alignment_sequence.length > longest_alignment
        longest_alignment = alignment.alignment_sequence.length
      end
      if AAsequence.all(:seq_id => alignment.seq_id, :contact_consensus.gte => 0.5).count > 0
        @cicp_contact_count += 1
      end
      alignment_array << alignment
    end
    for i in 0..longest_alignment+1
      @contact_consensus_array[i] = Array.new(@seq_contact_count, 0)
      @cicp_array[i] = Array.new(@seq_contact_count, 0)
    end
    #@contact_consensus_array = Array.new(longest_alignment, Array.new(@seq_contact_count,0))
    puts @contact_consensus_array.length
    puts "Into The Threads"
     thread_array=[]
      thread_num.times do |i|
        thread_array[i] = Thread.new{
          while alignment_array.length > 0 do
            alignment = alignment_array.pop
            sequence= alignment.sequence
            display_hash = Hash.new
            alignment_color_array = Array.new      
            cur_position = 0   
            orig_position = 0
            info_hash = {}
            seq_length = alignment.sequence.a_asequences.count
            contacts = alignment.sequence.a_asequences.all(:contact_consensus.gte =>0.6)
            disorders = alignment.sequence.a_asequences.all(:disorder_consensus.gte =>0.5)
            both = alignment.sequence.a_asequences.all(:disorder_consensus.gte =>0.5, :contact_consensus.gte => 0.6)
            info_hash[:abrev_name] = alignment.sequence.abrev_name
            info_hash[:cicp_range] = contacts.map{|aa| aa.original_position}.to_range
            info_hash[:cicp_counts] = contacts.count
            info_hash[:disorder_count] = disorders.count
            info_hash[:disorder_range] =  disorders.map{|aa| aa.original_position}.to_range
            info_hash[:disorder_percentage] = disorders.count.to_f/seq_length.to_f
            info_hash[:cicp_percentage] = contacts.count.to_f/seq_length.to_f
            info_hash[:both] = both.map{|aa| aa.original_position}.to_range
            info_hash[:both_count] = both.count
            info_hash[:both_percentage] = both.count.to_f/seq_length.to_f
            info_hash[:sequence_length] =seq_length
            @cicp_and_disorder_stats << info_hash
            AlignmentPosition.all(:alignment_id => alignment.align_id, 
                         :order => [:alignment_position_id.asc]).each do |position|  
             if position.position == cur_position
                amino_acid = sequence.a_asequences.first(:original_position=>orig_position) #AAsequence.first(:id => position.aasequence_id)
                unless amino_acid.nil?
                  alignment_color_array[cur_position] = residue_color(amino_acid.disorder_consensus, amino_acid.contact_consensus)
                  if @contact_consensus_array[cur_position][alignment.align_order].nil?
                    @contact_consensus_array[cur_position][alignment.align_order] = 0
                  end
                  if amino_acid.disorder_consensus >= 0.5
                    @contact_consensus_array[cur_position][alignment.align_order] = @contact_consensus_array[cur_position][alignment.align_order] + 1
                  end
                  if amino_acid.contact_consensus >= 0.5
                    @cicp_array[cur_position][alignment.align_order] = @cicp_array[cur_position][alignment.align_order] + 1
                  end
                else
                  puts "Amino Acid doesn't exits: #{sequence.abrev_name} | #{cur_position}:#{orig_position}" 
                  alignment_color_array[cur_position] = residue_color(0, 0)
                  @contact_consensus_array[cur_position][alignment.align_order] = 0
                  @cicp_array[cur_position][alignment.align_order] = 0
                end
             else
                while position.position > cur_position
                             alignment_color_array[cur_position] = "FFFFFF"
                             cur_position += 1
                end
                amino_acid = sequence.a_asequences.first(:original_position=>orig_position) #AAsequence.first(:id => position.aasequence_id)
                unless amino_acid.nil?
                  alignment_color_array[cur_position] = residue_color(amino_acid.disorder_consensus, amino_acid.contact_consensus)
                  if @contact_consensus_array[cur_position].nil?
                    puts "OH no " + alignment.sequence.abrev_name
                  end
                  if @contact_consensus_array[cur_position][alignment.align_order].nil?
                     @contact_consensus_array[cur_position][alignment.align_order] = 0
                  end
                  if amino_acid.disorder_consensus >= 0.5
                     @contact_consensus_array[cur_position][alignment.align_order] = @contact_consensus_array[cur_position][alignment.align_order] + 1
                  end
                  if amino_acid.contact_consensus >= 0.5
                   @cicp_array[cur_position][alignment.align_order] = @cicp_array[cur_position][alignment.align_order] + 1
                  end
                else
                    puts "Amino Acid doesn't exits: #{sequence.abrev_name} | #{cur_position}:#{orig_position}" 
                    alignment_color_array[cur_position] = residue_color(0, 0)
                    @contact_consensus_array[cur_position][alignment.align_order] = 0
                    @cicp_array[cur_position][alignment.align_order] = 0
                end
              end
              cur_position += 1
              orig_position +=1
              end                 
              puts display_hash["name"] = Sequence.first(:seq_id => alignment.seq_id).abrev_name 
              display_hash["alignment"] = alignment_color_array
              @display_array[alignment.align_order] = display_hash
            if @max_count < cur_position
                   @max_count = cur_position
            end
          end
        }
     end
     thread_array.map{|t| t.join}

     @contact_consensus_array = @contact_consensus_array.map{|a| a.inject(0){|sum,item| sum + item}}
     @cicp_array = @cicp_array.map{|a| a.inject(0){|sum,item| sum + item}}
     disorder_array = @contact_consensus_array.map{|dv| dv.to_f/@seq_contact_count}
     cicp_avgs = @cicp_array.map{|dv| dv.to_f/@cicp_contact_count}
     aa_array = Array(1..@contact_consensus_array.count)
     group_data = aa_array.zip(disorder_array, cicp_avgs,@conservation)
     require 'csv'
     @filename = "#{align.alignment_name}_display_data.csv"
     CSV.open("public/"+@filename, "w") do |csv|
       csv << ["Position","Disorder","CICP","Conservation"]
       group_data.each do |gd|
         csv << gd.map{|e| e.nil? ? 0 : e}
       end
     end
     @cicp_info=[]
     @cicp_info50=[]
     @cicp_info40=[]
     @cicp_info30=[]
     @dis_info = []
     @dis_info50 = []
     @dis_info40 = []
     @dis_info30 = []
     for i in 0..@cicp_array.length-1
      cp = @cicp_array[i].to_f/@cicp_contact_count
      d = @contact_consensus_array[i].to_f/@seq_contact_count
      if cp > 0
        @cicp_info << i
        if cp >= 0.5
          @cicp_info50 << i
        end
        if cp >= 0.4
          @cicp_info40 << i
        end
        if cp >= 0.3
          @cicp_info30 << i
        end
      end
      if d > 0
        @dis_info << i
        if d >= 0.5
          @dis_info50 << i
        end
        if d >= 0.4
          @dis_info40 << i
        end
        if d >= 0.3
          @dis_info30 << i
        end
      end
     end
     @cur_position = 0
     @tick_counter = 0
     @alignment_tick_array = Array.new
     while @cur_position <= @max_count
       @cur_position += 1
       @tick_counter += 1
       if @tick_counter != 25
         @alignment_tick_array << "FFFFFF"
       else
         @alignment_tick_array << "000000"
         @tick_counter = 0
       end
     end
     @display_hash = Hash.new
     @display_hash["name"] = ""
     @display_hash["alignment"] = @alignment_tick_array  
     @display_array << @display_hash
     if params[:aa_length].nil?
       @aa_length = 400
     else
       @aa_length = params[:aa_length].to_i
     end
     @ranges = (@max_count/@aa_length)

   end 


   def display_disorder_and_cicp_and_inter_annotated_alignment
      thread_num = 65
      @display_array = Array.new
      @max_count = 0
      @cicp_and_disorder_stats = []
      @contact_consensus_array = Array.new
      @cicp_array = Array.new
      @inter_consensus = Array.new
      @cicp_contact_count =0
      align = Alignment.get(params[:id])
      align2_ids = Alignment.get(params[:a2_id]).sequences.map{|s| s.seq_id}
      @seq_contact_count = Alignment.all(:alignment_name =>align.alignment_name).count
      longest_alignment = 0;
      alignment_array = []
      unless Conservation.first(:alignment_name => align.alignment_name)
        align.run_aacon
      end
      if Conservation.first(:alignment_name => align.alignment_name)
        @conservation = Conservation.first(:alignment_name => align.alignment_name).results_array
      else
        @conservation = []
      end
      @alignment_name = Alignment.get(params[:id]).alignment_name
      Alignment.all(:alignment_name => Alignment.get(params[:id]).alignment_name, 
                                  :order => [:align_order.asc]).each do |alignment|
       puts alignment.alignment_sequence.length
       if alignment.alignment_sequence.length > longest_alignment
         longest_alignment = alignment.alignment_sequence.length
       end
       if AAsequence.all(:seq_id => alignment.seq_id, :contact_consensus.gte => 0.5).count > 0
         @cicp_contact_count += 1
       end
       alignment_array << alignment
     end
     for i in 0..longest_alignment+1
       @contact_consensus_array[i] = Array.new(@seq_contact_count, 0)
       @cicp_array[i] = Array.new(@seq_contact_count, 0)
       @inter_consensus[i] = Array.new(@seq_contact_count, 0)
     end
     #@contact_consensus_array = Array.new(longest_alignment, Array.new(@seq_contact_count,0))
     puts @contact_consensus_array.length
     puts "Into The Threads"
      thread_array=[]
       thread_num.times do |i|
         thread_array[i] = Thread.new{
           while alignment_array.length > 0 do
             alignment = alignment_array.pop
             sequence= alignment.sequence
             display_hash = Hash.new
             alignment_color_array = Array.new      
             cur_position = 0   
             orig_position = 0
             info_hash = {}
             seq_length = alignment.sequence.a_asequences.count
             contacts = alignment.sequence.a_asequences.all(:contact_consensus.gte =>0.6)
             disorders = alignment.sequence.a_asequences.all(:disorder_consensus.gte =>0.5)
             both = alignment.sequence.a_asequences.all(:disorder_consensus.gte =>0.5, :contact_consensus.gte => 0.6)
             info_hash[:abrev_name] = alignment.sequence.abrev_name
             info_hash[:cicp_range] = contacts.map{|aa| aa.original_position}.to_range
             info_hash[:cicp_counts] = contacts.count
             info_hash[:disorder_count] = disorders.count
             info_hash[:disorder_range] =  disorders.map{|aa| aa.original_position}.to_range
             info_hash[:disorder_percentage] = disorders.count.to_f/seq_length.to_f
             info_hash[:cicp_percentage] = contacts.count.to_f/seq_length.to_f
             info_hash[:both] = both.map{|aa| aa.original_position}.to_range
             info_hash[:both_count] = both.count
             info_hash[:both_percentage] = both.count.to_f/seq_length.to_f
             info_hash[:sequence_length] =seq_length
             @cicp_and_disorder_stats << info_hash
             AlignmentPosition.all(:alignment_id => alignment.align_id, 
                          :order => [:alignment_position_id.asc]).each do |position|  
              if position.position == cur_position
                 amino_acid = sequence.a_asequences.first(:original_position=>orig_position) #AAsequence.first(:id => position.aasequence_id)
                 unless amino_acid.nil?
                   if InterCap.all(:aasequence1_id => amino_acid.id, :alignment1_id=>alignment.align_id, :seq2_id => align2_ids, :unique=>true, :fields=>[:aasequence1_id, :aasequence2_id]).count > 0  
                      intercap = InterCap.all(:aasequence1_id => amino_acid.id, :alignment1_id=>alignment.align_id, :seq2_id => align2_ids, :unique=>true, :fields=>[:aasequence1_id, :aasequence2_id]).count
                      @inter_consensus[cur_position][alignment.align_order] =1
                    elsif InterCap.all(:aasequence2_id => amino_acid.id, :seq1_id => align2_ids, :alignment2_id => alignment.align_id, :unique=>true, :fields=>[:aasequence1_id, :aasequence2_id]).count > 0  
                      intercap = InterCap.all(:aasequence2_id => amino_acid.id, :seq1_id => align2_ids, :alignment2_id => alignment.align_id, :unique=>true, :fields=>[:aasequence1_id, :aasequence2_id]).count
                      @inter_consensus[cur_position][alignment.align_order] =1
                    else
                      intercap = 0
                    end

                   alignment_color_array[cur_position] = inter_residue_color(amino_acid.disorder_consensus, amino_acid.contact_consensus,intercap)
                   if @contact_consensus_array[cur_position][alignment.align_order].nil?
                     @contact_consensus_array[cur_position][alignment.align_order] = 0
                   end
                   if amino_acid.disorder_consensus >= 0.5
                     @contact_consensus_array[cur_position][alignment.align_order] = @contact_consensus_array[cur_position][alignment.align_order] + 1
                   end
                   if amino_acid.contact_consensus >= 0.5
                     @cicp_array[cur_position][alignment.align_order] = @cicp_array[cur_position][alignment.align_order] + 1
                   end
                 else
                   puts "Amino Acid doesn't exits: #{sequence.abrev_name} | #{cur_position}:#{orig_position}" 
                   alignment_color_array[cur_position] = residue_color(0, 0)
                   @contact_consensus_array[cur_position][alignment.align_order] = 0
                   @cicp_array[cur_position][alignment.align_order] = 0
                 end
              else
                 while position.position > cur_position
                              alignment_color_array[cur_position] = "FFFFFF"
                              cur_position += 1
                 end
                 amino_acid = sequence.a_asequences.first(:original_position=>orig_position) #AAsequence.first(:id => position.aasequence_id)
                 unless amino_acid.nil?
                   if InterCap.all(:aasequence1_id => amino_acid.id, :alignment1_id=>alignment.align_id, :seq2_id => align2_ids, :unique=>true, :fields=>[:aasequence1_id, :aasequence2_id]).count > 0  
                       intercap = InterCap.all(:aasequence1_id => amino_acid.id,:alignment1_id=>alignment.align_id, :seq2_id => align2_ids,  :unique=>true, :fields=>[:aasequence1_id, :aasequence2_id]).count
                       @inter_consensus[cur_position][alignment.align_order] =1
                     elsif InterCap.all(:aasequence2_id => amino_acid.id, :seq1_id => align2_ids, :alignment2_id => alignment.align_id, :unique=>true, :fields=>[:aasequence1_id, :aasequence2_id]).count > 0  
                       intercap = InterCap.all(:aasequence2_id => amino_acid.id, :seq1_id => align2_ids, :alignment2_id => alignment.align_id, :unique=>true, :fields=>[:aasequence1_id, :aasequence2_id]).count
                       @inter_consensus[cur_position][alignment.align_order] =1
                     else
                       intercap = 0
                     end

                    alignment_color_array[cur_position] = inter_residue_color(amino_acid.disorder_consensus, amino_acid.contact_consensus,intercap)
                   if @contact_consensus_array[cur_position].nil?
                     puts "OH no " + alignment.sequence.abrev_name
                   end
                   if @contact_consensus_array[cur_position][alignment.align_order].nil?
                      @contact_consensus_array[cur_position][alignment.align_order] = 0
                   end
                   if amino_acid.disorder_consensus >= 0.5
                      @contact_consensus_array[cur_position][alignment.align_order] = @contact_consensus_array[cur_position][alignment.align_order] + 1
                   end
                   if amino_acid.contact_consensus >= 0.5
                    @cicp_array[cur_position][alignment.align_order] = @cicp_array[cur_position][alignment.align_order] + 1
                   end
                 else
                     puts "Amino Acid doesn't exits: #{sequence.abrev_name} | #{cur_position}:#{orig_position}" 
                     alignment_color_array[cur_position] = residue_color(0, 0)
                     @contact_consensus_array[cur_position][alignment.align_order] = 0
                     @cicp_array[cur_position][alignment.align_order] = 0
                 end
               end
               cur_position += 1
               orig_position +=1
               end                 
               puts display_hash["name"] = Sequence.first(:seq_id => alignment.seq_id).abrev_name 
               display_hash["alignment"] = alignment_color_array
               @display_array[alignment.align_order] = display_hash
             if @max_count < cur_position
                    @max_count = cur_position
             end
           end
         }
      end
      thread_array.map{|t| t.join}

      @contact_consensus_array = @contact_consensus_array.map{|a| a.inject(0){|sum,item| sum + item}}
      @inter_consensus = @inter_consensus.map{|a| a.inject(0){|sum,item| sum + item}}
      @cicp_array = @cicp_array.map{|a| a.inject(0){|sum,item| sum + item}}
      disorder_array = @contact_consensus_array.map{|dv| dv.to_f/@seq_contact_count}
      cicp_avgs = @cicp_array.map{|dv| dv.to_f/@cicp_contact_count}
      int_avgs =  @inter_consensus.map{|dv| dv.to_f/@cicp_contact_count}
      aa_array = Array(1..@contact_consensus_array.count)
      group_data = aa_array.zip(disorder_array, cicp_avgs,int_avgs,@conservation)
      require 'csv'
      @filename = "#{align.alignment_name}_inter_display_data.csv"
      CSV.open("public/"+@filename, "w") do |csv|
        csv << ["Position","Disorder","CICP","Inter", "Conservation"]
        group_data.each do |gd|
          csv << gd.map{|e| e.nil? ? 0 : e}
        end
      end
      @cicp_info=[]
      @cicp_info50=[]
      @cicp_info40=[]
      @cicp_info30=[]
      @dis_info = []
      @dis_info50 = []
      @dis_info40 = []
      @dis_info30 = []
      for i in 0..@cicp_array.length-1
       cp = @cicp_array[i].to_f/@cicp_contact_count
       d = @contact_consensus_array[i].to_f/@seq_contact_count
       if cp > 0
         @cicp_info << i
         if cp >= 0.5
           @cicp_info50 << i
         end
         if cp >= 0.4
           @cicp_info40 << i
         end
         if cp >= 0.3
           @cicp_info30 << i
         end
       end
       if d > 0
         @dis_info << i
         if d >= 0.5
           @dis_info50 << i
         end
         if d >= 0.4
           @dis_info40 << i
         end
         if d >= 0.3
           @dis_info30 << i
         end
       end
      end
      @cur_position = 0
      @tick_counter = 0
      @alignment_tick_array = Array.new
      while @cur_position <= @max_count
        @cur_position += 1
        @tick_counter += 1
        if @tick_counter != 25
          @alignment_tick_array << "FFFFFF"
        else
          @alignment_tick_array << "000000"
          @tick_counter = 0
        end
      end
      @display_hash = Hash.new
      @display_hash["name"] = ""
      @display_hash["alignment"] = @alignment_tick_array  
      @display_array << @display_hash
      if params[:aa_length].nil?
        @aa_length = 400
      else
        @aa_length = params[:aa_length].to_i
      end
      @ranges = (@max_count/@aa_length)

    end



   def plotCICP
     thread_num = 65
     @seq_contact_count = Alignment.all(:alignment_name => Alignment.get(params[:id]).alignment_name).count
     @dis_array = []
     alignment_array =[]
     @cicp_contact_count =0
     
     @alignment_name = Alignment.get(params[:id]).alignment_name
     Alignment.all(:alignment_name => Alignment.get(params[:id]).alignment_name, 
                                 :order => [:align_order.asc]).each do |alignment|
      puts alignment.alignment_sequence.length
      if AAsequence.all(:seq_id => alignment.seq_id, :contact_consensus.gte => 0.5).count > 0
        @cicp_contact_count += 1
      end
      alignment_array << alignment
    end
    puts "Into The Threads"
     thread_array=[]
      thread_num.times do |i|
        thread_array[i] = Thread.new{
          while alignment_array.length > 0 do
            alignment = alignment_array.pop
            sequence= alignment.sequence
            display_hash = Hash.new
            alignment_color_array = Array.new      
            cur_position = 0   
            orig_position = 0
            AlignmentPosition.all(:alignment_id => alignment.align_id, 
                         :order => [:alignment_position_id.asc]).each do |position|
             if position.position == cur_position
                amino_acid = sequence.a_asequences.first(:original_position=>orig_position)
                unless amino_acid.nil?
                  unless amino_acid.disorder_consensus.nil? || amino_acid.contact_consensus.nil?
                    @dis_array << "#{amino_acid.disorder_consensus},#{amino_acid.contact_consensus}"
                    #@dis_hash["#{amino_acid.disorder_consensus}"]["#{amino_acid.contact_consensus}"] = @dis_hash["#{amino_acid.disorder_consensus}"]["#{amino_acid.contact_consensus}"].nil? ? 1 : @dis_hash["#{amino_acid.disorder_consensus}"]["#{amino_acid.contact_consensus}"] + 1
                  end
                end
             else
                while position.position > cur_position
                             cur_position += 1
                end
                amino_acid = sequence.a_asequences.first(:original_position=>orig_position)
                unless amino_acid.nil?
                  unless amino_acid.disorder_consensus.nil? || amino_acid.contact_consensus.nil?
                    @dis_array << "#{amino_acid.disorder_consensus},#{amino_acid.contact_consensus}"
                    #@dis_hash["#{amino_acid.disorder_consensus}"]["#{amino_acid.contact_consensus}"] = @dis_hash["#{amino_acid.disorder_consensus}"]["#{amino_acid.contact_consensus}"].nil? ? 1 : @dis_hash["#{amino_acid.disorder_consensus}"]["#{amino_acid.contact_consensus}"] + 1
                  end
                end
             end
              cur_position += 1
              orig_position +=1          
            end
            puts Sequence.first(:seq_id => alignment.seq_id).abrev_name 
          end
        }
     end
     thread_array.map{|t| t.join}
   end 


   
   def display_cicp_annotated_alignment
     thread_num = 65
     @display_array = Array.new
     @max_count = 0
     @contact_consensus_array = Array.new
     @seq_contact_count = Alignment.all(:alignment_name => Alignment.get(params[:id]).alignment_name).count
     longest_alignment = 0;
     alignment_array = []
     @alignment_name = Alignment.get(params[:id]).alignment_name
     Alignment.all(:alignment_name => Alignment.get(params[:id]).alignment_name, 
                                 :order => [:align_order.asc]).each do |alignment|
      puts alignment.alignment_sequence.length
      if alignment.alignment_sequence.length > longest_alignment
        longest_alignment = alignment.alignment_sequence.length
      end
      alignment_array << alignment
    end
    for i in 0..longest_alignment
      @contact_consensus_array[i] = Array.new(@seq_contact_count, 0)
    end
    #@contact_consensus_array = Array.new(longest_alignment, Array.new(@seq_contact_count,0))
    puts @contact_consensus_array.length
    puts "Into The Threads"
     thread_array=[]
      thread_num.times do |i|
        thread_array[i] = Thread.new{
          while alignment_array.length > 0 do
            alignment = alignment_array.pop
            display_hash = Hash.new
            alignment_color_array = Array.new      
            cur_position = 0   
            AlignmentPosition.all(:alignment_id => alignment.align_id, 
                         :order => [:alignment_position_id.asc]).each do |position|
             if position.position == cur_position
                amino_acid = AAsequence.first(:id => position.aasequence_id)
                alignment_color_array[cur_position] = residue_color(amino_acid.disorder_consensus, amino_acid.contact_consensus)
                if @contact_consensus_array[cur_position][alignment.align_order].nil?
                  @contact_consensus_array[cur_position][alignment.align_order] = 0
                end
                if amino_acid.contact_consensus >= 0.5
                  @contact_consensus_array[cur_position][alignment.align_order] = @contact_consensus_array[cur_position][alignment.align_order] + 1
                end
             else
                while position.position > cur_position
                             alignment_color_array[cur_position] = "FFFFFF"
                             cur_position += 1
                end
                amino_acid = AAsequence.first(:id => position.aasequence_id)
                alignment_color_array[cur_position] = residue_color(amino_acid.disorder_consensus, amino_acid.contact_consensus)
                if @contact_consensus_array[cur_position][alignment.align_order].nil?
                   @contact_consensus_array[cur_position][alignment.align_order] = 0
                end
                if amino_acid.contact_consensus >= 0.5
                   @contact_consensus_array[cur_position][alignment.align_order] = @contact_consensus_array[cur_position][alignment.align_order] + 1
                end
              end
              cur_position += 1
              end                 
              puts display_hash["name"] = Sequence.first(:seq_id => alignment.seq_id).abrev_name 
              display_hash["alignment"] = alignment_color_array
              @display_array[alignment.align_order] = display_hash
            if @max_count < cur_position
                   @max_count = cur_position
            end
          end
        }
     end
     thread_array.map{|t| t.join}

     @contact_consensus_array = @contact_consensus_array.map{|a| a.inject(0){|sum,item| sum + item}}
     @cur_position = 0
     @tick_counter = 0
     @alignment_tick_array = Array.new
     while @cur_position <= @max_count
       @cur_position += 1
       @tick_counter += 1
       if @tick_counter != 25
         @alignment_tick_array << "FFFFFF"
       else
         @alignment_tick_array << "000000"
         @tick_counter = 0
       end
     end
     @display_hash = Hash.new
     @display_hash["name"] = ""
     @display_hash["alignment"] = @alignment_tick_array  
     @display_array << @display_hash
     if params[:aa_length].nil?
       @aa_length = 400
     else
       @aa_length = params[:aa_length].to_i
     end
     @ranges = (@max_count/@aa_length)

   end 



   
  
   def download_disorder_alignment
     thread_num = 65
     @display_array = Array.new
     @max_count = 0
     @contact_consensus_array = Array.new
     @seq_contact_count = Alignment.all(:alignment_name => Alignment.get(params[:id]).alignment_name).count
     longest_alignment = 0;
     alignment_array = []
     @alignment_name = Alignment.get(params[:id]).alignment_name
     Alignment.all(:alignment_name => Alignment.get(params[:id]).alignment_name, 
                                 :order => [:align_order.asc]).each do |alignment|
      puts alignment.alignment_sequence.length
      if alignment.alignment_sequence.length > longest_alignment
        longest_alignment = alignment.alignment_sequence.length
      end
      alignment_array << alignment
    end
    for i in 0..longest_alignment
      @contact_consensus_array[i] = Array.new(@seq_contact_count, 0)
    end
    #@contact_consensus_array = Array.new(longest_alignment, Array.new(@seq_contact_count,0))
    puts @contact_consensus_array.length
    puts "Into The Threads"
     thread_array=[]
      thread_num.times do |i|
        thread_array[i] = Thread.new{
          while alignment_array.length > 0 do
            alignment = alignment_array.pop
            display_hash = Hash.new
            alignment_color_array = Array.new      
            cur_position = 0   
            AlignmentPosition.all(:alignment_id => alignment.align_id, 
                         :order => [:alignment_position_id.asc]).each do |position|
             if position.position == cur_position
                amino_acid = AAsequence.first(:id => position.aasequence_id)
                alignment_color_array[cur_position] = residue_color(amino_acid.disorder_consensus, 0)
                if @contact_consensus_array[cur_position][alignment.align_order].nil?
                  @contact_consensus_array[cur_position][alignment.align_order] = 0
                end
                if amino_acid.disorder_consensus > 0.5
                  @contact_consensus_array[cur_position][alignment.align_order] = @contact_consensus_array[cur_position][alignment.align_order] + 1
                end
             else
                while position.position > cur_position
                             alignment_color_array[cur_position] = "FFFFFF"
                             cur_position += 1
                end
                amino_acid = AAsequence.first(:id => position.aasequence_id)
                alignment_color_array[cur_position] = residue_color(amino_acid.disorder_consensus, 0)
                if @contact_consensus_array[cur_position][alignment.align_order].nil?
                   @contact_consensus_array[cur_position][alignment.align_order] = 0
                end
                if amino_acid.disorder_consensus > 0.5
                   @contact_consensus_array[cur_position][alignment.align_order] = @contact_consensus_array[cur_position][alignment.align_order] + 1
                end
              end
              cur_position += 1
              end                 
              puts display_hash["name"] = Sequence.first(:seq_id => alignment.seq_id).abrev_name 
              display_hash["alignment"] = alignment_color_array
              @display_array[alignment.align_order] = display_hash
            if @max_count < cur_position
                   @max_count = cur_position
            end
          end
        }
      end
     thread_array.map{|t| t.join}

     @contact_consensus_array = @contact_consensus_array.map{|a| a.inject(0){|sum,item| sum + item}}
     @cur_position = 0
     @tick_counter = 0
     @alignment_tick_array = Array.new
     while @cur_position <= @max_count
       @cur_position += 1
       @tick_counter += 1
       if @tick_counter != 25
         @alignment_tick_array << "FFFFFF"
       else
         @alignment_tick_array << "000000"
         @tick_counter = 0
       end
     end
     @display_hash = Hash.new
     @display_hash["name"] = ""
     @display_hash["alignment"] = @alignment_tick_array  
     @display_array << @display_hash
     if params[:aa_length].nil?
       @aa_length = 400
     else
       @aa_length = params[:aa_length].to_i
     end
     @ranges = (@max_count/@aa_length)
     csv_string = CSV.generate do |csv|
       csv << ['position','disorder consensus'] 
       @contact_consensus_array.length.times do |i| 
         csv << [i.to_s ,  (@contact_consensus_array[i].to_f / @seq_contact_count).to_s]  
       end
     end
     send_data csv_string, 
     :type => 'text/csv; charset=iso-8859-1; header=present', 
     :disposition => "attachment; filename=#{@alignment_name}_disorder_consensus.csv"	
     flash[:notice] = "Export complete!" 
   end
     
  
  def display_annotated_alignment
    @display_array = Array.new
    @max_count = 0
    @contact_consensus_array = Array.new
    @seq_contact_count = 0
    Alignment.all(:alignment_name => Alignment.get(params[:id]).alignment_name, 
                                :order => [:align_order.asc]).each do |alignment|
      if AAsequence.all(:seq_id => alignment.seq_id, :contact_consensus.gt => 0.75).count > 0
        @seq_contact_count += 1
      end
      puts Sequence.first(:seq_id => alignment.seq_id).abrev_name + ":" + AAsequence.all(:seq_id => alignment.seq_id, :contact_consensus.gt => 0.75).count.to_s
    end
    
    Alignment.all(:alignment_name => Alignment.get(params[:id]).alignment_name, 
                                :order => [:align_order.asc]).each do |alignment|
    #@alignments = Alignment.all(:alignment_name => params[:aligment_name], 
                                #:order => [:align_order.asc]).each do |alignment|

       @display_hash = Hash.new
       @alignment_color_array = Array.new      
       @cur_position = 0   
       AlignmentPosition.all(:alignment_id => alignment.align_id, 
                    :order => [:alignment_position_id.asc]).each do |position|
        if position.position == @cur_position
           @amino_acid = AAsequence.first(:id => position.aasequence_id)
           @alignment_color_array[@cur_position] = residue_color(@amino_acid.disorder_consensus, @amino_acid.contact_consensus)
           if @contact_consensus_array[@cur_position].nil?
             @contact_consensus_array[@cur_position] = 0
           end
           if @amino_acid.contact_consensus > 0.75
             @contact_consensus_array[@cur_position] = @contact_consensus_array[@cur_position] + 1
           end
        else
           while position.position > @cur_position
                        @alignment_color_array[@cur_position] = "FFFFFF"
                        @cur_position += 1
           end
           @amino_acid = AAsequence.first(:id => position.aasequence_id)
           @alignment_color_array[@cur_position] = residue_color(@amino_acid.disorder_consensus, @amino_acid.contact_consensus)
           if @contact_consensus_array[@cur_position].nil?
              @contact_consensus_array[@cur_position] = 0
            end
           if @amino_acid.contact_consensus > 0.75
              @contact_consensus_array[@cur_position] = @contact_consensus_array[@cur_position] + 1
           end
         end
         @cur_position += 1
         end                 
         puts @display_hash["name"] = Sequence.first(:seq_id => alignment.seq_id).abrev_name 
         @display_hash["alignment"] = @alignment_color_array
         @display_array << @display_hash
       if @max_count < @cur_position
              @max_count = @cur_position
       end
    end
    @cur_position = 0
    @tick_counter = 0
    @alignment_tick_array = Array.new
    while @cur_position <= @max_count
      @cur_position += 1
      @tick_counter += 1
      if @tick_counter != 25
        @alignment_tick_array << "FFFFFF"
      else
        @alignment_tick_array << "000000"
        @tick_counter = 0
      end
    end
    @display_hash = Hash.new
    @display_hash["name"] = ""
    @display_hash["alignment"] = @alignment_tick_array  
    @display_array << @display_hash
    if params[:aa_length].nil?
      @aa_length = 400
    else
      @aa_length = params[:aa_length].to_i
    end
    @ranges = (@max_count/@aa_length)
    
  end
  
  def display_compensatory_annotated_alignment
    thread_num = 65
    @display_array = Array.new
    @max_count = 0
    @contact_consensus_array = Array.new
    @seq_contact_count = 0
    longest_alignment = 0
    orig_position = 0
    alignment_array = []
    alignment_count = Alignment.all(:alignment_name => Alignment.get(params[:id]).alignment_name).count
    @alignment_name = Alignment.get(params[:id]).alignment_name
    Alignment.all(:alignment_name => Alignment.get(params[:id]).alignment_name, 
                                :order => [:align_order.asc]).each do |alignment|
      if AAsequence.all(:seq_id => alignment.seq_id, :contact_consensus.gte => 0.5).count > 0
        @seq_contact_count += 1
      end
      puts Sequence.first(:seq_id => alignment.seq_id).abrev_name + ":" + AAsequence.all(:seq_id => alignment.seq_id, :contact_consensus.gte => 0.5).count.to_s
      if alignment.alignment_sequence.length > longest_alignment
        longest_alignment = alignment.alignment_sequence.length
      end
      alignment_array << alignment
    end
    #@contact_consensus_array = Array.new(longest_alignment,0)
    for i in 0..longest_alignment
      @contact_consensus_array[i] = Array.new(alignment_count, 0)
    end
    puts @contact_consensus_array.length
    puts "Into The Threads"
     thread_array=[]
      thread_num.times do |i|
        thread_array[i] = Thread.new{
          while alignment_array.length > 0 do
            alignment = alignment_array.pop
            sequence = alignment.sequence
            display_hash = Hash.new
            alignment_color_array = Array.new      
            cur_position = 0   
            AlignmentPosition.all(:alignment_id => alignment.align_id, 
                    :order => [:alignment_position_id.asc]).each do |position|
              if position.position == cur_position
                amino_acid = sequence.a_asequences.first(:original_position=>orig_position)#AAsequence.first(:id => position.aasequence_id)
                #cap_res = NewCap.first(:aasequence_id => amino_acid.AAsequence_id)
                cap_color =0
                # if amino_acid.contact_consensus >= 0.5 #@contact_consensus_array[@cur_position] > 1
                #                   cap_color = 1
                #                 # @contact_consensus_array[cur_position][alignment.align_order] = @contact_consensus_array[cur_position][alignment.align_order] + 1
                #                 end
                #                 @contact_consensus_array[cur_position][alignment.align_order] = @contact_consensus_array[cur_position][alignment.align_order] + amino_acid.contact_consensus
                #                 alignment_color_array[cur_position] = residue_color(amino_acid.disorder_consensus, cap_color)
                #                 
                if !amino_acid.nil? && amino_acid.contact_consensus >= 0.5 #@contact_consensus_array[@cur_position] > 1
                  cap_color = 1
                  @contact_consensus_array[cur_position][alignment.align_order] = @contact_consensus_array[cur_position][alignment.align_order] + amino_acid.contact_consensus
                 alignment_color_array[cur_position] = residue_color(amino_acid.disorder_consensus, cap_color)
                  # @contact_consensus_array[cur_position][alignment.align_order] = @contact_consensus_array[cur_position][alignment.align_order] + 1
                else
                  puts "Something Failed: " + position.alignment_position_id.to_s + " | position: " + position.position.to_s + " | aaseq_id: " + position.aasequence_id.to_s
                end
              else
               while position.position > cur_position
                alignment_color_array[cur_position] = "FFFFFF"
                cur_position += 1
               end
               amino_acid = sequence.a_asequences.first(:original_position=>orig_position)#AAsequence.first(:id => position.aasequence_id)
               #cap_res = NewCap.first(:aasequence_id => amino_acid.id)
                cap_color =0
                cap_color =0
                if !amino_acid.nil? && amino_acid.contact_consensus >= 0.5 #@contact_consensus_array[@cur_position] > 1
                  cap_color = 1
                  @contact_consensus_array[cur_position][alignment.align_order] = @contact_consensus_array[cur_position][alignment.align_order] + amino_acid.contact_consensus
                 alignment_color_array[cur_position] = residue_color(amino_acid.disorder_consensus, cap_color)
                  # @contact_consensus_array[cur_position][alignment.align_order] = @contact_consensus_array[cur_position][alignment.align_order] + 1
                else
                  puts "Something Failed: #{position.alignment_position_id.to_s} | position: #{position.position.to_s} | aaseq_id: #{position.aasequence_id.to_s}"
                end  
             end
           cur_position += 1
           orig_position += 1
         end                 
         puts display_hash["name"] = Sequence.first(:seq_id => alignment.seq_id).abrev_name 
         display_hash["alignment"] = alignment_color_array
         @display_array << display_hash
         if @max_count < cur_position
           @max_count = cur_position
         end
       end
       }
    end
    #@max_count=longest_alignment
    thread_array.map{|t| t.join}
    @contact_consensus_array = @contact_consensus_array.map{|a| a.inject(0){|sum,item| sum + item}}
    @cur_position = 0
    @tick_counter = 0
    @alignment_tick_array = Array.new
    while @cur_position <= @max_count
      @cur_position += 1
      @tick_counter += 1
      if @tick_counter != 25
        @alignment_tick_array << "FFFFFF"
      else
        @alignment_tick_array << "000000"
        @tick_counter = 0
      end
    end
    @display_hash = Hash.new
    @display_hash["name"] = ""
    @display_hash["alignment"] = @alignment_tick_array  
    @display_array << @display_hash
    if params[:aa_length].nil?
      @aa_length = 400
    else
      @aa_length = params[:aa_length].to_i
    end
    @ranges = (@max_count/@aa_length)
    if @seq_contact_count == 0
      @seq_contact_count = alignment_count
    end
  end


  def display_caps_annotated_alignment(thread_num=65)
    @display_array = Array.new
    @max_count = 0
    @contact_consensus_array = Array.new
    @seq_contact_count = 0
    longest_alignment = 0
    alignment_array = []
    alignment_count = Alignment.all(:alignment_name => Alignment.get(params[:id]).alignment_name).count
    @alignment_name = Alignment.get(params[:id]).alignment_name
    Alignment.all(:alignment_name => Alignment.get(params[:id]).alignment_name, 
                                :order => [:align_order.asc]).each do |alignment|
      if AAsequence.all(:seq_id => alignment.seq_id, :contact_consensus.gte => 0.5).count > 0
        @seq_contact_count += 1
      end
      puts Sequence.first(:seq_id => alignment.seq_id).abrev_name + ":" + AAsequence.all(:seq_id => alignment.seq_id, :contact_consensus.gte => 0.5).count.to_s
      if alignment.alignment_sequence.length > longest_alignment
        longest_alignment = alignment.alignment_sequence.length
      end
      alignment_array << alignment
    end
    #@contact_consensus_array = Array.new(longest_alignment,0)
    for i in 0..longest_alignment
      @contact_consensus_array[i] = Array.new(alignment_count, 0)
    end
    puts @contact_consensus_array.length
    puts "Into The Threads"
     thread_array=[]
      thread_num.times do |i|
        thread_array[i] = Thread.new{
          while alignment_array.length > 0 do
            alignment = alignment_array.pop
            display_hash = Hash.new
            alignment_color_array = Array.new      
            cur_position = 0   
            AlignmentPosition.all(:alignment_id => alignment.align_id, 
                    :order => [:alignment_position_id.asc]).each do |position|
              if position.position == @cur_position
                amino_acid = AAsequence.first(:id => position.aasequence_id)
                #cap_res = NewCap.first(:aasequence_id => amino_acid.AAsequence_id)
                cap_color =0
                if amino_acid.contact_consensus >= 0.5 #@contact_consensus_array[@cur_position] > 1
                  cap_color = 1
                end
                @contact_consensus_array[cur_position][alignment.align_order] = @contact_consensus_array[cur_position][alignment.align_order] + amino_acid.contact_consensus
                alignment_color_array[cur_position] = residue_color(amino_acid.disorder_consensus, cap_color)
              else
               while position.position > cur_position
                alignment_color_array[cur_position] = "FFFFFF"
                cur_position += 1
               end
               amino_acid = AAsequence.first(:id => position.aasequence_id)
               #cap_res = NewCap.first(:aasequence_id => amino_acid.id)
                cap_color =0
                cap_color =0
                if amino_acid.contact_consensus >= 0.5 #@contact_consensus_array[@cur_position] > 1
                  cap_color = 1
                end
                @contact_consensus_array[cur_position][alignment.align_order] = @contact_consensus_array[cur_position][alignment.align_order] + amino_acid.contact_consensus
               alignment_color_array[cur_position] = residue_color(amino_acid.disorder_consensus, cap_color)           
             end
           cur_position += 1
         end                 
         puts display_hash["name"] = Sequence.first(:seq_id => alignment.seq_id).abrev_name 
         display_hash["alignment"] = alignment_color_array
         @display_array << display_hash
         if @max_count < cur_position
           @max_count = cur_position
         end
       end
       }
    end
    #@max_count=longest_alignment
    thread_array.map{|t| t.join}
    @contact_consensus_array = @contact_consensus_array.map{|a| a.inject(0){|sum,item| sum + item}}
    @cur_position = 0
    @tick_counter = 0
    @alignment_tick_array = Array.new
    while @cur_position <= @max_count
      @cur_position += 1
      @tick_counter += 1
      if @tick_counter != 25
        @alignment_tick_array << "FFFFFF"
      else
        @alignment_tick_array << "000000"
        @tick_counter = 0
      end
    end
    @display_hash = Hash.new
    @display_hash["name"] = ""
    @display_hash["alignment"] = @alignment_tick_array  
    @display_array << @display_hash
    if params[:aa_length].nil?
      @aa_length = 400
    else
      @aa_length = params[:aa_length].to_i
    end
    @ranges = (@max_count/@aa_length)
    if @seq_contact_count == 0
      @seq_contact_count = alignment_count
    end
  end

   
  def display_xdet_annotated_alignment
     thread_num = 65
     @display_array = Array.new
     @max_count = 0
     @contact_consensus_array = Array.new
     @seq_contact_count = 0
     longest_alignment = 0
     alignment_array = []
     alignment_count = Alignment.all(:alignment_name => Alignment.get(params[:id]).alignment_name).count
     @alignment_name = Alignment.get(params[:id]).alignment_name
     Alignment.all(:alignment_name => Alignment.get(params[:id]).alignment_name, 
                                 :order => [:align_order.asc]).each do |alignment|
       if AAsequence.all(:seq_id => alignment.seq_id, :contact_consensus.gte => 0.5).count > 0
         @seq_contact_count += 1
       end
       puts Sequence.first(:seq_id => alignment.seq_id).abrev_name + ":" + AAsequence.all(:seq_id => alignment.seq_id, :contact_consensus.gt => 0.75).count.to_s
       if alignment.alignment_sequence.length > longest_alignment
         longest_alignment = alignment.alignment_sequence.length
       end
       alignment_array << alignment
     end
     #@contact_consensus_array = Array.new(longest_alignment,0)
     for i in 0..longest_alignment
       @contact_consensus_array[i] = Array.new(alignment_count, 0)
     end
     puts @contact_consensus_array.length
     puts "Into The Threads"
      thread_array=[]
       thread_num.times do |i|
         thread_array[i] = Thread.new{
           while alignment_array.length > 0 do
             alignment = alignment_array.pop
             display_hash = Hash.new
             alignment_color_array = Array.new      
             cur_position = 0   
             AlignmentPosition.all(:alignment_id => alignment.align_id, 
                     :order => [:alignment_position_id.asc]).each do |position|
               if position.position == @cur_position
                 amino_acid = AAsequence.first(:id => position.aasequence_id)
                 cap_res = Caps.first(:aasequence_id => amino_acid.AAsequence_id)
                 cap_color =0
                 if !amino_acid.xdet.nil?  #&& (amino_acid.xdet.correlation > 0.0 || amino_acid.xdet.correlation == -2) #@contact_consensus_array[@cur_position] > 1
                   cap_color = 1
                 # @contact_consensus_array[cur_position][alignment.align_order] = @contact_consensus_array[cur_position][alignment.align_order] + 1
                 end
                 @contact_consensus_array[cur_position][alignment.align_order] = @contact_consensus_array[cur_position][alignment.align_order] + cap_color
                 alignment_color_array[cur_position] = residue_color(amino_acid.disorder_consensus, cap_color)
               else
                while position.position > cur_position
                 alignment_color_array[cur_position] = "FFFFFF"
                 cur_position += 1
                end
                amino_acid = AAsequence.first(:id => position.aasequence_id)
                cap_res = Caps.first(:aasequence_id => amino_acid.AAsequence_id)
                 cap_color =0
                 cap_color =0
                 if !amino_acid.xdet.nil? # && (amino_acid.xdet.correlation > 0.0 || amino_acid.xdet.correlation == -2)#@contact_consensus_array[@cur_position] > 1
                   cap_color = 1
                   # @contact_consensus_array[cur_position][alignment.align_order] = @contact_consensus_array[cur_position][alignment.align_order] + 1
                 end
                 @contact_consensus_array[cur_position][alignment.align_order] = @contact_consensus_array[cur_position][alignment.align_order] + cap_color
                alignment_color_array[cur_position] = residue_color(amino_acid.disorder_consensus, cap_color)           
              end
            cur_position += 1
          end                 
          puts display_hash["name"] = Sequence.first(:seq_id => alignment.seq_id).abrev_name 
          display_hash["alignment"] = alignment_color_array
          @display_array << display_hash
          if @max_count < cur_position
            @max_count = cur_position
          end
        end
        }
     end
     #@max_count=longest_alignment
     thread_array.map{|t| t.join}
     @contact_consensus_array = @contact_consensus_array.map{|a| a.inject(0){|sum,item| sum + item}}
     @cur_position = 0
     @tick_counter = 0
     @alignment_tick_array = Array.new
     while @cur_position <= @max_count
       @cur_position += 1
       @tick_counter += 1
       if @tick_counter != 25
         @alignment_tick_array << "FFFFFF"
       else
         @alignment_tick_array << "000000"
         @tick_counter = 0
       end
     end
     @display_hash = Hash.new
     @display_hash["name"] = ""
     @display_hash["alignment"] = @alignment_tick_array  
     @display_array << @display_hash
     if params[:aa_length].nil?
       @aa_length = 400
     else
       @aa_length = params[:aa_length].to_i
     end
     @ranges = (@max_count/@aa_length)
     if @seq_contact_count == 0
       @seq_contact_count = alignment_count
     end
   end
   
   def display_caps_annotated_alignment
      thread_num = 65
      @display_array = Array.new
      @max_count = 0
      @contact_consensus_array = Array.new
      @seq_contact_count = 0
      longest_alignment = 0
      alignment_array = []
      alignment_count = Alignment.all(:alignment_name => Alignment.get(params[:id]).alignment_name).count
      @alignment_name = Alignment.get(params[:id]).alignment_name
      Alignment.all(:alignment_name => Alignment.get(params[:id]).alignment_name, 
                                  :order => [:align_order.asc]).each do |alignment|
        if AAsequence.all(:seq_id => alignment.seq_id, :contact_consensus.gte => 0.5).count > 0
          @seq_contact_count += 1
        end
        puts Sequence.first(:seq_id => alignment.seq_id).abrev_name + ":" + AAsequence.all(:seq_id => alignment.seq_id, :contact_consensus.gt => 0.75).count.to_s
        if alignment.alignment_sequence.length > longest_alignment
          longest_alignment = alignment.alignment_sequence.length
        end
        alignment_array << alignment
      end
      #@contact_consensus_array = Array.new(longest_alignment,0)
      for i in 0..longest_alignment
        @contact_consensus_array[i] = Array.new(alignment_count, 0)
      end
      puts @contact_consensus_array.length
      puts "Into The Threads"
       thread_array=[]
        thread_num.times do |i|
          thread_array[i] = Thread.new{
            while alignment_array.length > 0 do
              alignment = alignment_array.pop
              display_hash = Hash.new
              alignment_color_array = Array.new      
              cur_position = 0   
              AlignmentPosition.all(:alignment_id => alignment.align_id, 
                      :order => [:alignment_position_id.asc]).each do |position|
                if position.position == @cur_position
                  amino_acid = AAsequence.first(:id => position.aasequence_id)
                  cap_res = Caps.first(:aasequence_id => amino_acid.AAsequence_id)
                  cap_color =0
                  if !NewCap.first(:seq_id=> amino_acid.seq_id, :position_one => amino_acid.original_position).nil? || !NewCap.first(:seq_id=> amino_acid.seq_id, :position_two => amino_acid.original_position).nil?
                    cap_color = 1
                  end

                  @contact_consensus_array[cur_position][alignment.align_order] = @contact_consensus_array[cur_position][alignment.align_order] + cap_color
                  alignment_color_array[cur_position] = residue_color(amino_acid.disorder_consensus, cap_color)
                else
                 while position.position > cur_position
                  alignment_color_array[cur_position] = "FFFFFF"
                  cur_position += 1
                 end
                 amino_acid = AAsequence.first(:id => position.aasequence_id)
                 cap_res = Caps.first(:aasequence_id => amino_acid.AAsequence_id)
                  cap_color =0
                  cap_color =0
                  if !NewCap.first(:seq_id=> amino_acid.seq_id, :position_one => amino_acid.original_position).nil? || !NewCap.first(:seq_id=> amino_acid.seq_id, :position_two => amino_acid.original_position).nil?
                    cap_color = 1
                  end
                  @contact_consensus_array[cur_position][alignment.align_order] = @contact_consensus_array[cur_position][alignment.align_order] + cap_color
                 alignment_color_array[cur_position] = residue_color(amino_acid.disorder_consensus, cap_color)           
               end
             cur_position += 1
           end                 
           puts display_hash["name"] = Sequence.first(:seq_id => alignment.seq_id).abrev_name 
           display_hash["alignment"] = alignment_color_array
           @display_array << display_hash
           if @max_count < cur_position
             @max_count = cur_position
           end
         end
         }
      end
      #@max_count=longest_alignment
      thread_array.map{|t| t.join}
      @contact_consensus_array = @contact_consensus_array.map{|a| a.inject(0){|sum,item| sum + item}}
      @cur_position = 0
      @tick_counter = 0
      @alignment_tick_array = Array.new
      while @cur_position <= @max_count
        @cur_position += 1
        @tick_counter += 1
        if @tick_counter != 25
          @alignment_tick_array << "FFFFFF"
        else
          @alignment_tick_array << "000000"
          @tick_counter = 0
        end
      end
      @display_hash = Hash.new
      @display_hash["name"] = ""
      @display_hash["alignment"] = @alignment_tick_array  
      @display_array << @display_hash
      if params[:aa_length].nil?
        @aa_length = 400
      else
        @aa_length = params[:aa_length].to_i
      end
      @ranges = (@max_count/@aa_length)
      if @seq_contact_count == 0
        @seq_contact_count = alignment_count
      end
    end


    
  
  def caps_report(thread_num=70)
    @caps_array = []
    alignment = Alignment.get(params[:id])
    @alignment_name = alignment.alignment_name
    seq_array = alignment.sequences.to_a
    thread_array=[]
    thread_num.times do |i|
      thread_array[i] = Thread.new{
        while seq_array.length > 0 do
          s = seq_array.pop
          temp_hash = {}
          temp_hash[:name] = s.abrev_name
          temp_hash[:id] = s.seq_id
          temp_hash[:caps_count] = s.new_caps_count
          temp_hash["20"] = s.new_caps_count_gt_twenty
          temp_hash["50"] = s.new_caps_count_gt_50
          temp_hash["100"] = s.new_caps_count_gt_100
          temp_hash["200"] = s.new_caps_count_gt_200
          temp_hash["300"] = s.new_caps_count_gt_300
          temp_hash["400"] = s.new_caps_count_gt_400
          temp_hash["500"] = s.new_caps_count_gt_500
          temp_hash["600"] = s.new_caps_count_gt_600
          temp_hash["700"] = s.new_caps_count_gt_700
          temp_hash["800"] = s.new_caps_count_gt_800
          temp_hash["900"] = s.new_caps_count_gt_900
          temp_hash["1000"] = s.new_caps_count_gt_1000
          temp_hash["1100"] = s.new_caps_count_gt_1100
          temp_hash["1200"] = s.new_caps_count_gt_1200
          temp_hash["1300"] = s.new_caps_count_gt_1300
          temp_hash["1400"] = s.new_caps_count_gt_1400
          temp_hash["1500"] = s.new_caps_count_gt_1500
          temp_hash["1600"] = s.new_caps_count_gt_1600
          temp_hash["1700"] = s.new_caps_count_gt_1700
          temp_hash["1800"] = s.new_caps_count_gt_1800
          temp_hash["1900"] = s.new_caps_count_gt_1900
          temp_hash["2000"] = s.new_caps_count_gt_2000
          temp_hash["2100"] = s.new_caps_count_gt_2100
          temp_hash["2200"] = s.new_caps_count_gt_2200
          temp_hash["2300"] = s.new_caps_count_gt_2300
          temp_hash["2400"] = s.new_caps_count_gt_2400
          temp_hash["2500"] = s.new_caps_count_gt_2500
          temp_hash[:pids] = PercentIdentity.all(:seq1_id => s.seq_id, :percent_id.gt => 19,:percent_id.lt => 90, :alignment_name => alignment.alignment_name).count
          @caps_array << temp_hash
        end
      }
    end
    #@max_count=longest_alignment
    thread_array.map{|t| t.join}
  end
  
  def inter_residue_color(dis_avg, con_avg, inter)
    if inter > 0
      if con_avg > 0.6
        if dis_avg >= 0.5
         @color = "FA1D2F"  #raspberry red, dis, con and inter
        else
         @color =  "E066FF" #medium orchid1, con and inter
        end
      elsif dis_avg >= 0.5
        @color ="FFCC00" #orange both dis and inter
      else
        @color = "0000FF" # blue "FF9966" #peach, inter only
      end
    elsif con_avg > 0.75
      if dis_avg >= 0.5
       @color = "00FF00"  #3394560  3CF3C0 green , dis and con
      else
       @color =  "66FFCC" #6750156  66FFCC, light blue con
      end
    else  #color for disorder only
      if dis_avg >= 0.5 && dis_avg < 0.6  #yellow for disorder only
       @color = "FFFF00" #16776960  FFFF00
      elsif dis_avg >= 0.6 && dis_avg < 0.7
       @color =  "FFFF00" #16763904  FFCC00
      elsif dis_avg >= 0.7 && dis_avg < 0.8
       @color =  "FFFF00" #16750848  FF9900
      elsif dis_avg >= 0.8 && dis_avg < 0.9
       @color = "FFFF00"  #16737792  FF6600
      elsif dis_avg >= 0.9
       @color =  "FFFF00" #16711680  FF0000
      else
       @color = "999999"
      end
    end
    return @color
  end
end
