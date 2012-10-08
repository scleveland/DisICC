class XdetsController < ApplicationController
  # GET /xdets
  # GET /xdets.xml
  def index
    @xdets = Xdet.all

    respond_to do |format|
      format.html # index.html.erb
      format.xml  { render :xml => @xdets }
    end
  end

  # GET /xdets/1
  # GET /xdets/1.xml
  def show
    @xdet = Xdet.find(params[:id])

    respond_to do |format|
      format.html # show.html.erb
      format.xml  { render :xml => @xdet }
    end
  end

  # GET /xdets/new
  # GET /xdets/new.xml
  def new
    @xdet = Xdet.new

    respond_to do |format|
      format.html # new.html.erb
      format.xml  { render :xml => @xdet }
    end
  end

  # GET /xdets/1/edit
  def edit
    @xdet = Xdet.find(params[:id])
  end

  # POST /xdets
  # POST /xdets.xml
  def create
    @xdet = Xdet.new(params[:xdet])

    respond_to do |format|
      if @xdet.save
        format.html { redirect_to(@xdet, :notice => 'Xdet was successfully created.') }
        format.xml  { render :xml => @xdet, :status => :created, :location => @xdet }
      else
        format.html { render :action => "new" }
        format.xml  { render :xml => @xdet.errors, :status => :unprocessable_entity }
      end
    end
  end

  # PUT /xdets/1
  # PUT /xdets/1.xml
  def update
    @xdet = Xdet.find(params[:id])

    respond_to do |format|
      if @xdet.update_attributes(params[:xdet])
        format.html { redirect_to(@xdet, :notice => 'Xdet was successfully updated.') }
        format.xml  { head :ok }
      else
        format.html { render :action => "edit" }
        format.xml  { render :xml => @xdet.errors, :status => :unprocessable_entity }
      end
    end
  end

  # DELETE /xdets/1
  # DELETE /xdets/1.xml
  def destroy
    @xdet = Xdet.find(params[:id])
    @xdet.destroy

    respond_to do |format|
      format.html { redirect_to(xdets_url) }
      format.xml  { head :ok }
    end
  end
end
