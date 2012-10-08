class ConseqsController < ApplicationController
  # GET /conseqs
  # GET /conseqs.xml
  def index
    @conseqs = Conseq.all

    respond_to do |format|
      format.html # index.html.erb
      format.xml  { render :xml => @conseqs }
    end
  end

  # GET /conseqs/1
  # GET /conseqs/1.xml
  def show
    @conseq = Conseq.find(params[:id])

    respond_to do |format|
      format.html # show.html.erb
      format.xml  { render :xml => @conseq }
    end
  end

  # GET /conseqs/new
  # GET /conseqs/new.xml
  def new
    @conseq = Conseq.new

    respond_to do |format|
      format.html # new.html.erb
      format.xml  { render :xml => @conseq }
    end
  end

  # GET /conseqs/1/edit
  def edit
    @conseq = Conseq.find(params[:id])
  end

  # POST /conseqs
  # POST /conseqs.xml
  def create
    @conseq = Conseq.new(params[:conseq])

    respond_to do |format|
      if @conseq.save
        format.html { redirect_to(@conseq, :notice => 'Conseq was successfully created.') }
        format.xml  { render :xml => @conseq, :status => :created, :location => @conseq }
      else
        format.html { render :action => "new" }
        format.xml  { render :xml => @conseq.errors, :status => :unprocessable_entity }
      end
    end
  end

  # PUT /conseqs/1
  # PUT /conseqs/1.xml
  def update
    @conseq = Conseq.find(params[:id])

    respond_to do |format|
      if @conseq.update_attributes(params[:conseq])
        format.html { redirect_to(@conseq, :notice => 'Conseq was successfully updated.') }
        format.xml  { head :ok }
      else
        format.html { render :action => "edit" }
        format.xml  { render :xml => @conseq.errors, :status => :unprocessable_entity }
      end
    end
  end

  # DELETE /conseqs/1
  # DELETE /conseqs/1.xml
  def destroy
    @conseq = Conseq.find(params[:id])
    @conseq.destroy

    respond_to do |format|
      format.html { redirect_to(conseqs_url) }
      format.xml  { head :ok }
    end
  end
  
  def import_conseq_data(dir, seq_type)
    Dir.foreach(dir) do |file|
       if file.count(".") == 0
         myfile = File.new(dir+'/'+file, "r")
         line_num = 1
         abrev_name = file.split("_")[0]
         puts abrev_name+'*********************************'
         seq = Sequence.first(:abrev_name => abrev_name,:seq_type => seq_type)
         while (line = myfile.gets)
           if line_num > 13 && !line.rstrip.empty?
             results = line.split("\t")
             #POS	 SEQ	 SCORE	COLOR	B/E	FUNCTION	MSA DATA	RESIDUE VARIETY
             puts results[0].lstrip+':'+results[1].lstrip+':'+results[2].lstrip+':'+results[3].lstrip+':'+results[4].lstrip+':'+results[5].lstrip+':'+results[6].lstrip+':'+results[7].lstrip
             if !AAsequence.first(:seq_id => seq.seq_id, :original_position => results[0].lstrip.to_i).nil?
               conseq = Conseq.first_or_create(:seq_id => seq.seq_id,
                                   :aasequence_id => AAsequence.first(:seq_id => seq.seq_id, :original_position => results[0].lstrip.to_i).AAsequence_id,
                                   :score => results[2].lstrip.to_f,
                                   :color => results[3].lstrip.to_i,
                                   :state => results[4].lstrip,
                                   :function => results[5].lstrip,
                                   :msa_data => results[6].lstrip,
                                   :residue_variety => results[7].lstrip)
             end
           end
           line_num +=1
         end
       end
    end
  end
end
