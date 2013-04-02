class InterCapsController < ApplicationController
  # GET /inter_caps
  # GET /inter_caps.xml
  def index
    sql = "Select Distinct seq1_id, seq2_id From inter_caps"
    results = repository.adapter.select(sql)
    @inter_caps = []
    results.each do |res|
      @inter_caps << InterCap.first(:seq1_id => res.seq1_id, :seq2_id => res.seq2_id)
    end

    respond_to do |format|
      format.html # index.html.erb
      format.xml  { render :xml => @inter_caps }
    end
  end

  # GET /inter_caps/1
  # GET /inter_caps/1.xml
  def show
    @inter_cap = InterCap.get(params[:id])
    @display_array = []
    @max_count = 0
    sequence_array =[@inter_cap.seq1,@inter_cap.seq2]
    puts "Into The Threads"
    thread_array=[]
    2.times do |i|
      thread_array[i] = Thread.new{
        while sequence_array.length > 0 do
          sequence= sequence_array.pop
          alignment_color_array = []
          contact_consensus_array= []
          disorder_consensus_array= []
          inter_array= []
          interflag_array =[]
          interflag_range=[]
          interflagfalse_range=[]
          display_hash = Hash.new
          alignment_color_array = Array.new      
          cur_position = 0   
          orig_position = 0
          info_hash = {}
          seq_length = sequence.a_asequences.count
          contacts = sequence.a_asequences.all(:contact_consensus.gte =>0.6)
          disorders = sequence.a_asequences.all(:disorder_consensus.gte =>0.5)
          both = sequence.a_asequences.all(:disorder_consensus.gte =>0.5, :contact_consensus.gte => 0.6)
          info_hash[:abrev_name] = sequence.abrev_name
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
          
          sequence.a_asequences.all(:order => [:original_position]).each do |amino_acid| 
           if InterCap.all(:aasequence1_id => amino_acid.id, :unique=>true, :fields=>[:aasequence1_id, :aasequence2_id]).count > 0  
             intercap = InterCap.all(:aasequence1_id => amino_acid.id, :unique=>true, :fields=>[:aasequence1_id, :aasequence2_id]).count
             interflag = 1
             interflag_range << amino_acid.original_position
           elsif InterCap.all(:aasequence2_id => amino_acid.id, :unique=>true, :fields=>[:aasequence1_id, :aasequence2_id]).count > 0  
             intercap = InterCap.all(:aasequence2_id => amino_acid.id, :unique=>true, :fields=>[:aasequence1_id, :aasequence2_id]).count
             interflag = 1
             interflag_range << amino_acid.original_position
           else
             intercap = 0
             inter_flag = 0
             interflagfalse_range << amino_acid.original_position
           end
           alignment_color_array[amino_acid.original_position] = residue_color(amino_acid.disorder_consensus, amino_acid.contact_consensus, intercap)
           contact_consensus_array[amino_acid.original_position]=amino_acid.contact_consensus
           disorder_consensus_array[amino_acid.original_position]=amino_acid.disorder_consensus
           inter_array[amino_acid.original_position]= intercap
           interflag_array[amino_acid.original_position]= interflag
          end 
          puts display_hash["name"] = sequence.abrev_name 
          display_hash["alignment"] = alignment_color_array
          display_hash["contacts"] = contact_consensus_array
          display_hash["disorders"]  = disorder_consensus_array
          display_hash["inter"] = inter_array
          display_hash["interflag"] = interflag_array
          display_hash["info"] = info_hash
          display_hash["inter_range"] = interflag_range
          display_hash["interfalse_range"] = interflagfalse_range
          @display_array << display_hash
        end
      }
    end
    thread_array.map{|t| t.join}
    # 
    # @contact_consensus_array = @contact_consensus_array.map{|a| a.inject(0){|sum,item| sum + item}}
    # @cicp_array = @cicp_array.map{|a| a.inject(0){|sum,item| sum + item}}
    # disorder_array = @contact_consensus_array.map{|dv| dv.to_f/@seq_contact_count}
    # cicp_avgs = @cicp_array.map{|dv| dv.to_f/@cicp_contact_count}
    # aa_array = Array(1..@contact_consensus_array.count)
    # group_data = aa_array.zip(disorder_array, cicp_avgs,@conservation)
    # require 'csv'
    # @filename = "#{align.alignment_name}_display_data.csv"
    # CSV.open("public/"+@filename, "w") do |csv|
    #  csv << ["Position","Disorder","CICP","Conservation"]
    #  group_data.each do |gd|
    #    csv << gd.map{|e| e.nil? ? 0 : e}
    #  end
    # end
    #debugger
    cur_position = 0
    tick_counter = 0
    alignment_tick_array = Array.new
    @max_count = @inter_cap.seq1.a_asequences.count > @inter_cap.seq2.a_asequences.count ? @inter_cap.seq1.a_asequences.count : @inter_cap.seq2.a_asequences.count
    
    while cur_position <= @max_count
     cur_position += 1
     tick_counter += 1
     if tick_counter != 25
       alignment_tick_array << "FFFFFF"
     else
       alignment_tick_array << "000000"
       tick_counter = 0
     end
    end
    display_hash = Hash.new
    display_hash["name"] = ""
    display_hash["alignment"] = alignment_tick_array  
    @display_array << display_hash
    if params[:aa_length].nil?
     @aa_length = 400
    else
     @aa_length = params[:aa_length].to_i
    end
    @ranges = (@max_count/@aa_length)

  
    respond_to do |format|
      format.html # show.html.erb
      format.xml  { render :xml => @inter_cap }
    end
  end

  # GET /inter_caps/new
  # GET /inter_caps/new.xml
  def new
    @inter_cap = InterCap.new

    respond_to do |format|
      format.html # new.html.erb
      format.xml  { render :xml => @inter_cap }
    end
  end

  # GET /inter_caps/1/edit
  def edit
    @inter_cap = InterCap.find(params[:id])
  end

  # POST /inter_caps
  # POST /inter_caps.xml
  def create
    @inter_cap = InterCap.new(params[:inter_cap])

    respond_to do |format|
      if @inter_cap.save
        format.html { redirect_to(@inter_cap, :notice => 'A asequence was successfully created.') }
        format.xml  { render :xml => @inter_cap, :status => :created, :location => @inter_cap }
      else
        format.html { render :action => "new" }
        format.xml  { render :xml => @inter_cap.errors, :status => :unprocessable_entity }
      end
    end
  end

  # PUT /inter_caps/1
  # PUT /inter_caps/1.xml
  def update
    @inter_cap = InterCap.find(params[:id])

    respond_to do |format|
      if @inter_cap.update_attributes(params[:inter_cap])
        format.html { redirect_to(@inter_cap, :notice => 'A asequence was successfully updated.') }
        format.xml  { head :ok }
      else
        format.html { render :action => "edit" }
        format.xml  { render :xml => @inter_cap.errors, :status => :unprocessable_entity }
      end
    end
  end

  # DELETE /inter_caps/1
  # DELETE /inter_caps/1.xml
  def destroy
    @inter_cap = InterCap.find(params[:id])
    @inter_cap.destroy

    respond_to do |format|
      format.html { redirect_to(inter_caps_url) }
      format.xml  { head :ok }
    end
  end
  
  private
  
  def residue_color(dis_avg, con_avg, inter)
    if inter > 0
      if con_avg > 0.6
        if dis_avg >= 0.5
         @color = "FA1D2F"  #raspberry red
        else
         @color =  "E066FF" #medium orchid1
        end
      elsif dis_avg >= 0.5
        @color ="2E0854" #indigo
      else
        @color = "FF00CC" #rose
      end
    elsif con_avg > 0.75
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
  
end
