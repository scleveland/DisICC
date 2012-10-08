class AasequencesController < ApplicationController
  # GET /aasequences
  # GET /aasequences.xml
  def index
    @aasequences = Aasequences.find(:all)

    respond_to do |format|
      format.html # index.html.erb
      format.xml  { render :xml => @aasequences }
    end
  end

  # GET /aasequences/1
  # GET /aasequences/1.xml
  def show
    @aasequences = Aasequences.find(params[:id])

    respond_to do |format|
      format.html # show.html.erb
      format.xml  { render :xml => @aasequences }
    end
  end

  # GET /aasequences/new
  # GET /aasequences/new.xml
  def new
    @aasequences = Aasequences.new

    respond_to do |format|
      format.html # new.html.erb
      format.xml  { render :xml => @aasequences }
    end
  end

  # GET /aasequences/1/edit
  def edit
    @aasequences = Aasequences.find(params[:id])
  end

  # POST /aasequences
  # POST /aasequences.xml
  def create
    @aasequences = Aasequences.new(params[:aasequences])

    respond_to do |format|
      if @aasequences.save
        flash[:notice] = 'Aasequences was successfully created.'
        format.html { redirect_to(@aasequences) }
        format.xml  { render :xml => @aasequences, :status => :created, :location => @aasequences }
      else
        format.html { render :action => "new" }
        format.xml  { render :xml => @aasequences.errors, :status => :unprocessable_entity }
      end
    end
  end

  # PUT /aasequences/1
  # PUT /aasequences/1.xml
  def update
    @aasequences = Aasequences.find(params[:id])

    respond_to do |format|
      if @aasequences.update_attributes(params[:aasequences])
        flash[:notice] = 'Aasequences was successfully updated.'
        format.html { redirect_to(@aasequences) }
        format.xml  { head :ok }
      else
        format.html { render :action => "edit" }
        format.xml  { render :xml => @aasequences.errors, :status => :unprocessable_entity }
      end
    end
  end

  # DELETE /aasequences/1
  # DELETE /aasequences/1.xml
  def destroy
    @aasequences = Aasequences.find(params[:id])
    @aasequences.destroy

    respond_to do |format|
      format.html { redirect_to(aasequences_url) }
      format.xml  { head :ok }
    end
  end
end
