class DisordersController < ApplicationController
  # GET /disorders
  # GET /disorders.xml
  def index
    @disorders = Disorder.all

    respond_to do |format|
      format.html # index.html.erb
      format.xml  { render :xml => @disorders }
    end
  end

  # GET /disorders/1
  # GET /disorders/1.xml
  def show
    @disorder = Disorder.find(params[:id])

    respond_to do |format|
      format.html # show.html.erb
      format.xml  { render :xml => @disorder }
    end
  end

  # GET /disorders/new
  # GET /disorders/new.xml
  def new
    @disorder = Disorder.new

    respond_to do |format|
      format.html # new.html.erb
      format.xml  { render :xml => @disorder }
    end
  end

  # GET /disorders/1/edit
  def edit
    @disorder = Disorder.find(params[:id])
  end

  # POST /disorders
  # POST /disorders.xml
  def create
    @disorder = Disorder.new(params[:disorder])

    respond_to do |format|
      if @disorder.save
        format.html { redirect_to(@disorder, :notice => 'Disorder was successfully created.') }
        format.xml  { render :xml => @disorder, :status => :created, :location => @disorder }
      else
        format.html { render :action => "new" }
        format.xml  { render :xml => @disorder.errors, :status => :unprocessable_entity }
      end
    end
  end

  # PUT /disorders/1
  # PUT /disorders/1.xml
  def update
    @disorder = Disorder.find(params[:id])

    respond_to do |format|
      if @disorder.update_attributes(params[:disorder])
        format.html { redirect_to(@disorder, :notice => 'Disorder was successfully updated.') }
        format.xml  { head :ok }
      else
        format.html { render :action => "edit" }
        format.xml  { render :xml => @disorder.errors, :status => :unprocessable_entity }
      end
    end
  end

  # DELETE /disorders/1
  # DELETE /disorders/1.xml
  def destroy
    @disorder = Disorder.find(params[:id])
    @disorder.destroy

    respond_to do |format|
      format.html { redirect_to(disorders_url) }
      format.xml  { head :ok }
    end
  end
end
