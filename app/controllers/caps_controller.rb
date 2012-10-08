class CapsController < ApplicationController
  # GET /caps
  # GET /caps.xml
  def index
    @caps = Cap.all

    respond_to do |format|
      format.html # index.html.erb
      format.xml  { render :xml => @caps }
    end
  end

  # GET /caps/1
  # GET /caps/1.xml
  def show
    @cap = Cap.find(params[:id])

    respond_to do |format|
      format.html # show.html.erb
      format.xml  { render :xml => @cap }
    end
  end

  # GET /caps/new
  # GET /caps/new.xml
  def new
    @cap = Cap.new

    respond_to do |format|
      format.html # new.html.erb
      format.xml  { render :xml => @cap }
    end
  end

  # GET /caps/1/edit
  def edit
    @cap = Cap.find(params[:id])
  end

  # POST /caps
  # POST /caps.xml
  def create
    @cap = Cap.new(params[:cap])

    respond_to do |format|
      if @cap.save
        format.html { redirect_to(@cap, :notice => 'Cap was successfully created.') }
        format.xml  { render :xml => @cap, :status => :created, :location => @cap }
      else
        format.html { render :action => "new" }
        format.xml  { render :xml => @cap.errors, :status => :unprocessable_entity }
      end
    end
  end

  # PUT /caps/1
  # PUT /caps/1.xml
  def update
    @cap = Cap.find(params[:id])

    respond_to do |format|
      if @cap.update_attributes(params[:cap])
        format.html { redirect_to(@cap, :notice => 'Cap was successfully updated.') }
        format.xml  { head :ok }
      else
        format.html { render :action => "edit" }
        format.xml  { render :xml => @cap.errors, :status => :unprocessable_entity }
      end
    end
  end

  # DELETE /caps/1
  # DELETE /caps/1.xml
  def destroy
    @cap = Cap.find(params[:id])
    @cap.destroy

    respond_to do |format|
      format.html { redirect_to(caps_url) }
      format.xml  { head :ok }
    end
  end
  
  def import_caps_data(dir, seq_type)
    Dir.foreach(dir) do |file|
       if file.split(".")[1] == "out"
         puts file
         myfile = File.new(dir + '/' + file, "r")
         line_num = 1
         abrev_name = file.split(".")[0].to_s
         puts abrev_name + '*********************************'
         seq = Sequence.first(:abrev_name => abrev_name,:seq_type => seq_type)
         start_line = 999999999
         while (line = myfile.gets)
           if line_num > start_line && !line.rstrip.empty?
             if line.count("Groups") ==1
               start_line = 999999999
             else
               results = line.split("\t")
               puts results.length
               puts line
               #Posicion AA1	Posicion AA2	Mean D1		Mean D2		Correlation
               puts results[1]# + ':' + results[2]#+':'+results[2]#.lstrip+':'+results[4].lstrip#+':'+results[6].lstrip #+':'+results[8].lstrip+':'+results[6].lstrip+':'+results[7].lstrip
               # if !AAsequence.first(:seq_id => seq.seq_id, :original_position => results[0].lstrip.to_i).nil?
               #   conseq = Conseq.first_or_create(:seq_id => seq.seq_id,
               #                       :aasequence_id => AAsequence.first(:seq_id => seq.seq_id, :original_position => results[0].lstrip.to_i).AAsequence_id,
               #                       :score => results[2].lstrip.to_f,
               #                       :color => results[3].lstrip.to_i,
               #                       :state => results[4].lstrip,
               #                       :function => results[5].lstrip,
               #                       :msa_data => results[6].lstrip,
               #                       :residue_variety => results[7].lstrip)
               # end
             end
           end
           line_num +=1
           if line.count("Posicion") > 1
             start_line = line_num +1
           end
         end
       end
    end
  end
end
