class SequencesController < ApplicationController
  # GET /sequences
  # GET /sequences.xml
  def index
    @sequences = Sequence.all#current_user.sequences#Sequence.all

    respond_to do |format|
      format.html # index.html.erb
      format.xml  { render :xml => @sequences }
    end
  end

  # GET /sequences/1
  # GET /sequences/1.xml
  def show
    @sequence = Sequence.get(params[:id])#current_user.sequences.get(params[:id])#Sequence.get(params[:id])

    respond_to do |format|
      format.html # show.html.erb
      format.xml  { render :xml => @sequence }
    end
  end

  # GET /sequences/new
  # GET /sequences/new.xml
  def new
    @sequence = Sequence.new
    @seq_types = {:N=>"N", :L=>"L", :P => "P", "P Cut" => "P Cut"}
    @sequences = Hash.new
    Sequence.all(:fields => [:seq_name], :unique => true, :order => [:seq_name.asc]).each{|s| @sequences = @sequences.merge({s.seq_name => s.seq_name})}
    respond_to do |format|
      format.html # new.html.erb
      format.xml  { render :xml => @sequence }
    end
  end

  # GET /sequences/1/edit
  def edit
    @sequence = Sequence.get(params[:id])#current_user.sequences.get(params[:id])#Sequence.get(params[:id])
  end

  # POST /sequences
  # POST /sequences.xml
  def create
    @sequence = Sequence.new(params[:sequence])
    @sequence.owner = current_user.id
    respond_to do |format|
      if @sequence.save
        begin
          @sequence.run_and_store_disorder()
        rescue Exception => e  
          puts @sequence.abrev_name + " Diorder Generation Failed"
          puts e.message
        end
        begin
          @sequence.calculate_disorder_consensus()
        rescue Exception => e  
          puts @sequence.abrev_name + " Diorder Consensus Calculation Failed"
          puts e.message
        end
        format.html { redirect_to(@sequence, :notice => 'Sequence was successfully created.') }
        format.xml  { render :xml => @sequence, :status => :created, :location => @sequence }
      else
        format.html { render :action => "new" }
        format.xml  { render :xml => @sequence.errors, :status => :unprocessable_entity }
      end
    end
  end

  # PUT /sequences/1
  # PUT /sequences/1.xml
  def update
    @sequence = Sequence.first(:seq_id=>params[:id].to_i)#current_user.sequences.first(:seq_id=>params[:id].to_i)# Sequence.get(params[:id])

    respond_to do |format|
      if @sequence.update_attributes(params[:sequence])
        format.html { redirect_to(@sequence, :notice => 'Sequence was successfully updated.') }
        format.xml  { head :ok }
      else
        format.html { render :action => "edit" }
        format.xml  { render :xml => @sequence.errors, :status => :unprocessable_entity }
      end
    end
  end

  # DELETE /sequences/1
  # DELETE /sequences/1.xml
  def destroy
    
    @sequence = Sequence.first(:seq_id=>params[:id].to_i) #current_user.sequences.first(:seq_id=>params[:id].to_i)#Sequence.get(params[:id])
    debugger
    @sequence.destroy!

    respond_to do |format|
      format.html { redirect_to(sequences_url) }
      format.xml  { head :ok }
    end
  end
  
  def disorder_consensus
    @sequence = Sequence.get(params[:id])#current_user.sequences.get(params[:id])
    respond_to do |format|
      format.html 
    end
  end
  
  def download_disorder
    @sequence = current_user.sequences.get(params[:id])
    csv_string = CSV.generate do |csv|
      csv << ['position','amino acid','disorder consensus'] + @sequence.disorders.all(:order=>[:disorder_type]).map{|d| d.disorder_type}
      AAsequence.all(:seq_id=>@sequence.seq_id, :order=>[:original_position]).each do |aa|
        csv << [(aa.original_position+1).to_s,aa.amino_acid,aa.disorder_consensus] + @sequence.disorders.all(:order=>[:disorder_type]).map{|d| d.disorder_values.first(:a_asequence_id=>aa.id).dvalue}
      end
    end

    filename = @sequence.seq_name + ".csv"
    send_data(csv_string,
    :type => 'text/csv; charset=utf-8; header=present',
    :filename => filename)
  end
  
  def run_disorder
    @sequence = Sequence.first(:seq_id=>params[:id].to_i)#current_user.sequences.first(:seq_id=>params[:id].to_i)
    @sequence.run_and_store_disorder()
    @sequence.calculate_disorder_consensus()
    redirect_to(sequences_url)
  end
end
