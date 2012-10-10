class SequencesController < ApplicationController
  # GET /sequences
  # GET /sequences.xml
  def index
    @sequences = current_user.sequences#Sequence.all

    respond_to do |format|
      format.html # index.html.erb
      format.xml  { render :xml => @sequences }
    end
  end

  # GET /sequences/1
  # GET /sequences/1.xml
  def show
    @sequence = current_user.sequences.get(params[:id])#Sequence.get(params[:id])

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
    @sequence = current_user.sequences.get(params[:id])#Sequence.get(params[:id])
  end

  # POST /sequences
  # POST /sequences.xml
  def create
    @sequence = Sequence.new(params[:sequence])
    @sequence.owner = current_user.id
    respond_to do |format|
      if @sequence.save
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
    @sequence = current_user.sequences.get(params[:id])# Sequence.get(params[:id])

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
    @sequence = current_user.sequences.get(params[:id])#Sequence.get(params[:id])
    @sequence.destroy

    respond_to do |format|
      format.html { redirect_to(sequences_url) }
      format.xml  { head :ok }
    end
  end
end
