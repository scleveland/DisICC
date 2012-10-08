class AlignmentsController < ApplicationController
  # GET /alignments
  # GET /alignments.xml
  def index
    @alignments = Alignment.all

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
    @alignment = Alignment.find(params[:id])
    @alignment.destroy

    respond_to do |format|
      format.html { redirect_to(alignments_url) }
      format.xml  { head :ok }
    end
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
  
  def pre_process_fasta_file
    @name = Time.now.to_s + params[:datafile].original_filename
    directory = "temp_data"
    @new_file = File.join(directory,@name)
    @alignment_name = params[:alignment_name]
    @seq_type = params[:seq_type]
  
      sequences = Sequence.all(:seq_type => @seq_type)
    
    @seq_options = Hash.new
    sequences.map{|k| @seq_options = @seq_options.merge({k.seq_name => k.seq_id.to_s})}
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
    @new_file = File.join(directory,params[:datafile_name])
    logger.debug "BEFORE*************************"
    begin
      file = File.new(@new_file, "r")
      logger.debug "After FILE OPEN"
      puts "After FILE OPEN"
      counter = 0
      order_count = 0
      abrev_name = ""
      alignment_sequence = ""
      fasta_hash = Hash.new
      logger.debug "BEFORE FASTA HASH*****"
      puts "BEFORE FASTA HASH*****"
      logger.debug "BEFORE@1232132131****"
      range = params[:seq_num].to_i - 1
      logger.debug "BEFORE FASTA RANGE***********"
      puts "BEFORE FASTA RANGE***********"
      (0..range).each do |i|
        logger.debug "BLAH"
        fasta_hash = fasta_hash.merge(params["fasta_name"+i.to_s].strip => params["seq"+i.to_s])
      end
      logger.debug {fasta_hash}
      logger.debug "AFTER FASTA HASH*******"
      while (line = file.gets)
        if line.count(">") > 0 && counter > 0
          #save the current sequence to an alignment
          logger.debug {fasta_hash[abrev_name]}
          puts fasta_hash[abrev_name]
          @sequence = Sequence.get(fasta_hash[abrev_name.strip])
          logger.debug "After"
          puts "After"  
          logger.debug "OHNO NONONONONO" 
          @alignment = Alignment.new(:seq_id => @sequence.seq_id,
                           :alignment_name => params[:alignment_name],
                           :align_order => order_count,
                           :alignment_sequence => alignment_sequence,
                           :fasta_title => abrev_name)
           logger.debug "SHAZZZZZAAAAAAAM"                 
          logger.debug { @alignment.valid?}
          if !@alignment.valid?
            puts @alignment.errors.inspect()
          end
          logger.debug "VALID"
          logger.debug { @alignment.errors.inspect }
          @alignment.save
          alignment_to_positions(@alignment)              
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
       puts fasta_hash[abrev_name]
        @sequence = Sequence.get(fasta_hash[abrev_name.strip])
        logger.debug "After"
        puts "After"  
        logger.debug "OHNO NONONONONO" 
        @alignment = Alignment.new(:seq_id => @sequence.seq_id,
                         :alignment_name => params[:alignment_name],
                         :align_order => order_count,
                         :alignment_sequence => alignment_sequence,
                         :fasta_title => abrev_name)
         logger.debug "SHAZZZZZAAAAAAAM"                 
        logger.debug { @alignment.valid?}
        if !@alignment.valid?
          puts @alignment.errors.inspect()
        end
        logger.debug "VALID"
        logger.debug { @alignment.errors.inspect }
        @alignment.save
        alignment_to_positions(@alignment)              
        #this is the sequene label
        abrev_name = line.gsub(">", "")

        order_count += 1
        alignment_sequence = ""
      file.close
      rescue => err
          puts "Exception: #{err}"
          err
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
          @alignment = Alignment.new(:seq_id => @sequence.seq_id,
                           :alignment_name => params[:alignment_name],
                           :align_order => order_count,
                           :alignment_sequence => alignment_sequence)
          @alignment.valid?
          logger.debug "VALID"
          logger.debug { @alignment.errors.inspect }
          @alignment.save
          alignment_to_positions(@alignment)              
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
        if !Caps.first(:seq_id=> aaseq.seq_id, :position_one => aaseq.original_position).nil?
          count +=1
        elsif !Caps.first(:seq_id=> aaseq.seq_id, :position_two => aaseq.original_position).nil?
          count +=1
        end
        aaseq.contact_positive_consensus = count /4
        aaseq.save
      end  
    end                          
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
     @display_array = Array.new
     @max_count = 0
     @contact_consensus_array = Array.new
     @seq_contact_count = 0
     Alignment.all(:alignment_name => Alignment.get(params[:id]).alignment_name, 
                                 :order => [:align_order.asc]).each do |alignment|
       #if AAsequence.all(:seq_id => alignment.seq_id, :contact_consensus.gt => 0.75).count > 0
         @seq_contact_count += 1
       #end
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
            @amino_acid = AAsequence.first(:AAsequence_id => position.aasequence_id)
            @alignment_color_array[@cur_position] = residue_color(@amino_acid.disorder_consensus, 0)
            if @contact_consensus_array[@cur_position].nil?
              @contact_consensus_array[@cur_position] = 0
            end
            if @amino_acid.disorder_consensus > 0.5
              @contact_consensus_array[@cur_position] = @contact_consensus_array[@cur_position] + 1
            end
         else
            while position.position > @cur_position
                         @alignment_color_array[@cur_position] = "FFFFFF"
                         @cur_position += 1
            end
            @amino_acid = AAsequence.first(:AAsequence_id => position.aasequence_id)
            @alignment_color_array[@cur_position] = residue_color(@amino_acid.disorder_consensus, 0)
            if @contact_consensus_array[@cur_position].nil?
               @contact_consensus_array[@cur_position] = 0
             end
            if @amino_acid.disorder_consensus > 0.5
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
           @amino_acid = AAsequence.first(:AAsequence_id => position.aasequence_id)
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
           @amino_acid = AAsequence.first(:AAsequence_id => position.aasequence_id)
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
           @amino_acid = AAsequence.first(:AAsequence_id => position.aasequence_id)
           @cap_res = Caps.first(:aasequence_id => @amino_acid.AAsequence_id)
           cap_color =0
           if !@cap_res.nil?
             cap_color = 1
           end
           @alignment_color_array[@cur_position] = residue_color(@amino_acid.disorder_consensus, cap_color)
           if @contact_consensus_array[@cur_position].nil?
             @contact_consensus_array[@cur_position] = 0
           end
           if !@cap_res.nil?
             @contact_consensus_array[@cur_position] = @contact_consensus_array[@cur_position] + 1
           end
        else
           while position.position > @cur_position
                        @alignment_color_array[@cur_position] = "FFFFFF"
                        @cur_position += 1
           end
           @amino_acid = AAsequence.first(:AAsequence_id => position.aasequence_id)
           @cap_res = Caps.first(:aasequence_id => @amino_acid.AAsequence_id)
            cap_color =0
            if !@cap_res.nil?
              cap_color = 1
            end
           @alignment_color_array[@cur_position] = residue_color(@amino_acid.disorder_consensus, cap_color)
           if @contact_consensus_array[@cur_position].nil?
              @contact_consensus_array[@cur_position] = 0
            end
           if !@cap_res.nil?
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
  
end
