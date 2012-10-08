class AAsequencesController < ApplicationController
  # GET /a_asequences
  # GET /a_asequences.xml
  def index
    @a_asequences = AAsequence.all

    respond_to do |format|
      format.html # index.html.erb
      format.xml  { render :xml => @a_asequences }
    end
  end

  # GET /a_asequences/1
  # GET /a_asequences/1.xml
  def show
    @a_asequence = AAsequence.find(params[:id])

    respond_to do |format|
      format.html # show.html.erb
      format.xml  { render :xml => @a_asequence }
    end
  end

  # GET /a_asequences/new
  # GET /a_asequences/new.xml
  def new
    @a_asequence = AAsequence.new

    respond_to do |format|
      format.html # new.html.erb
      format.xml  { render :xml => @a_asequence }
    end
  end

  # GET /a_asequences/1/edit
  def edit
    @a_asequence = AAsequence.find(params[:id])
  end

  # POST /a_asequences
  # POST /a_asequences.xml
  def create
    @a_asequence = AAsequence.new(params[:a_asequence])

    respond_to do |format|
      if @a_asequence.save
        format.html { redirect_to(@a_asequence, :notice => 'A asequence was successfully created.') }
        format.xml  { render :xml => @a_asequence, :status => :created, :location => @a_asequence }
      else
        format.html { render :action => "new" }
        format.xml  { render :xml => @a_asequence.errors, :status => :unprocessable_entity }
      end
    end
  end

  # PUT /a_asequences/1
  # PUT /a_asequences/1.xml
  def update
    @a_asequence = AAsequence.find(params[:id])

    respond_to do |format|
      if @a_asequence.update_attributes(params[:a_asequence])
        format.html { redirect_to(@a_asequence, :notice => 'A asequence was successfully updated.') }
        format.xml  { head :ok }
      else
        format.html { render :action => "edit" }
        format.xml  { render :xml => @a_asequence.errors, :status => :unprocessable_entity }
      end
    end
  end

  # DELETE /a_asequences/1
  # DELETE /a_asequences/1.xml
  def destroy
    @a_asequence = AAsequence.find(params[:id])
    @a_asequence.destroy

    respond_to do |format|
      format.html { redirect_to(a_asequences_url) }
      format.xml  { head :ok }
    end
  end
end
